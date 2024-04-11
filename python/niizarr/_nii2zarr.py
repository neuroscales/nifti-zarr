import sys
import json
import math
import zarr.hierarchy
import zarr.storage
import numpy as np
import numcodecs
import argparse
from skimage.transform import pyramid_gaussian, pyramid_laplacian
from ome_zarr.writer import write_multiscale
from nibabel import Nifti1Image, load
from ._header import (
    UNITS, DTYPES, INTENTS, INTENTS_P, SLICEORDERS, XFORMS,
    HEADERTYPE1, HEADERTYPE2,
)


# If fsspec available, use fsspec
try:
    import fsspec
    open = fsspec.open
except (ImportError, ModuleNotFoundError):
    fsspec = None


SYS_BYTEORDER = '<' if sys.byteorder == 'little' else '>'


def nii2json(header):
    """
    Convert a nifti header to JSON

    Parameters
    ----------
    header : np.array[HEADERTYPE1 or HEADERTYPE2]
        Nifti header in binary form

    Returns
    -------
    header : dict
        Nifti header in JSON form
    """
    # use nz1/nz2 instead of n+1/n+2
    header = header.copy()
    magic = header["magic"].tobytes().decode()
    magic = magic[:1] + 'z' + magic[2:]
    header['magic'] = magic

    ndim = header["dim"][0].item()
    intent_code = INTENTS[header["intent_code"].item()]
    jsonheader = {
        "magic": header["magic"].tobytes().decode(),
        "dim": header["dim"][1:1+ndim].tolist(),
        "pixdim": header["pixdim"][1:1+ndim].tolist(),
        "units": {
            "space": UNITS[(header["xyzt_units"] & 0x07).item()],
            "time": UNITS[(header["xyzt_units"] & 0x38).item()],
        },
        "datatype": DTYPES[header["datatype"].item()],
        "dim_info": {
            "freq": (header["dim_info"] & 0x03).item(),
            "phase": ((header["dim_info"] >> 2) & 0x03).item(),
            "slice": ((header["dim_info"] >> 4) & 0x03).item(),
        },
        "intent": {
            "code": intent_code,
            "name": header["intent_name"].tobytes().decode(),
            "p": header["intent_p"][:INTENTS_P[intent_code]].tolist(),
        },
        "scl": {
            "slope": header["scl_slope"].item(),
            "inter": header["scl_inter"].item(),
        },
        "slice": {
            "code": SLICEORDERS[header["slice_code"].item()],
            "start": header["slice_start"].item(),
            "end": header["slice_end"].item(),
            "duration": header["slice_duration"].item(),
        },
        "cal": {
            "min": header["cal_min"].item(),
            "max": header["cal_max"].item(),
        },
        "toffset": header["toffset"].item(),
        "descrip": header["descrip"].tobytes().decode(),
        "aux_file": header["aux_file"].tobytes().decode(),
        "qform": {
            "intent": XFORMS[header["qform_code"].item()],
            "quatern": header["quatern"].tolist(),
            "offset": header["qoffset"].tolist(),
            "fac": header["pixdim"][0].item(),
        },
        "sform": {
            "intent": XFORMS[header["sform_code"].item()],
            "affine": header["sform"].tolist(),
        },
    }
    if not math.isfinite(jsonheader["scl"]["slope"]):
        jsonheader["scl"]["slope"] = 0.0
    if not math.isfinite(jsonheader["scl"]["inter"]):
        jsonheader["scl"]["inter"] = 0.0

    # Fix data type
    byteorder = header['sizeof_hdr'].dtype.byteorder
    if byteorder == '=':
        byteorder = SYS_BYTEORDER
    if isinstance(jsonheader["datatype"], tuple):
        jsonheader["datatype"] = [
            [field, '|' + dtype] for field, dtype in jsonheader["datatype"]
        ]
    elif jsonheader["datatype"].endswith('1'):
        jsonheader["datatype"] = '|' + jsonheader["datatype"]
    else:
        jsonheader["datatype"] = byteorder + jsonheader["datatype"]

    # Check that the dictionary is serializable
    json.dumps(jsonheader)

    return jsonheader


def _make_compressor(name, **prm):
    if not isinstance(name, str):
        return name
    name = name.lower()
    if name == 'blosc':
        Compressor = numcodecs.Blosc
    elif name == 'zlib':
        Compressor = numcodecs.Zlib
    else:
        raise ValueError('Unknown compressor', name)
    return Compressor(**prm)


def _make_pyramid3d(
        data3d, nb_levels, pyramid_fn=pyramid_gaussian, label=False,
        no_pyramid_axis=None,
):
    no_pyramid_axis = {
        'x': 0,
        'y': 1,
        'z': 2,
    }.get(no_pyramid_axis, no_pyramid_axis)
    if isinstance(no_pyramid_axis, str):
        no_pyramid_axis = int(no_pyramid_axis)

    batch, nxyz = data3d.shape[:-3], data3d.shape[-3:]
    data3d = data3d.reshape((-1, *nxyz))
    max_layer = nb_levels - 1 if nb_levels > 0 else -1

    def pyramid_values(x):
        return pyramid_fn(x, max_layer, 2, preserve_range=True,
                          channel_axis=no_pyramid_axis)

    def pyramid_labels(x):
        yield x
        labels = np.unique(x)
        pyrmaxprob = list(pyramid_values(x == labels[0]))[1:]
        pyramid = [
            np.zeros_like(level, dtype=x.dtype) for level in pyrmaxprob
        ]
        for label in labels[1:]:
            pyrprob = list(pyramid_values(x == label))[1:]
            for (value, prob, maxprob) in zip(pyramid, pyrprob, pyrmaxprob):
                mask = prob > maxprob
                value[mask] = label
                maxprob[mask] = prob[mask]

        for level in pyramid:
            yield level

    pyramid = pyramid_labels if label else pyramid_values

    for level in zip(*map(pyramid, data3d)):
        yield np.stack(level).reshape(batch + level[0].shape)


def nii2zarr(inp, out, *,
             chunk=64,
             chunk_channel=1,
             chunk_time=1,
             nb_levels=-1,
             method='gaussian',
             label=None,
             no_time=False,
             no_pyramid_axis=None,
             fill_value=None,
             compressor='blosc',
             compressor_options={}):
    """
    Convert a nifti file to nifti-zarr

    Parameters
    ----------
    inp : nib.Nifti1Image or file_like
        Input nifti image
    out : zarr.Store or zarr.Group or path
        Output zarr object

    Other Parameters
    ----------------
    chunk : [list of [tuple of]] int
        Chunk size of the spatial dimensions, per x/y/z, per level.
        * The inner tuple allows different chunk sizes to be used along
          each dimension.
        * The outer list allows different chunk sizes to be used at
          different pyramid levels.
    chunk_channel : int
        Chunk size of the channel dimension. If 0, combine all channels
        in a single chunk.
    chunk_time : int
        Chunk size of the time dimension. If 0, combine all timepoints
        in a single chunk.
    nb_levels : int
        Number of levels in the pyramid. Default: all possible levels.
    method : {'gaussian', 'laplacian'}
        Method used to compute the pyramid.
    label : bool
        Is this is a label volume?  If `None`, guess from intent code.
    no_time : bool
        If True, there is no time dimension so the 4th dimension
        (if it exists) should be interpreted as the channel dimensions.
    no_pyramid_axis : int or {'x', 'y', 'z'}
        Axis that should not be downsampled. If None, downsample
        across all three dimensions.
    fill_value : number
        Value to use for missing tiles
    compressor : {'blosc', 'zlib'}
        Compression to use
    compressor_options : dict
        Compressor options
    """
    # Open nifti image with nibabel
    if not isinstance(inp, Nifti1Image):
        if hasattr(inp, 'read'):
            inp = Nifti1Image.from_stream(inp)
        else:
            inp = load(inp)

    # Open zarr group
    if not isinstance(out, zarr.hierarchy.Group):
        if not isinstance(out, zarr.storage.Store):
            if fsspec:
                out = zarr.storage.FSStore(out)
            else:
                out = zarr.storage.DirectoryStore(out)
        out = zarr.group(store=out, overwrite=True)

    if no_time and len(inp.shape) > 3:
        inp = Nifti1Image(inp.dataobj[:, :, :, None], inp.affine, inp.header)

    # Parse nifti header
    v = int(inp.header.structarr['magic'].tobytes().decode()[2])
    header = np.frombuffer(inp.header.structarr.tobytes(), count=1,
                           dtype=HEADERTYPE1 if v == 1 else HEADERTYPE2)[0]
    if header['magic'].decode() not in ('ni1', 'n+1', 'ni2', 'n+2'):
        swappedheader = header.newbyteorder()
    else:
        swappedheader = header
    jsonheader = nii2json(swappedheader)

    data = np.asarray(inp.dataobj)

    if fill_value:
        if np.issubdtype(data.dtype, np.complexfloating):
            fill_value = complex(fill_value)
        elif np.issubdtype(data.dtype, np.floating):
            fill_value = float(fill_value)
        elif np.issubdtype(data.dtype, np.integer):
            fill_value = int(fill_value)
        elif np.issubdtype(data.dtype, np.bool_):
            fill_value = bool(fill_value)

    # Fix array shape
    if data.ndim == 5:
        nbatch = 2
        perm = [3, 4, 2, 1, 0]
        axes = ['t', 'c', 'z', 'y', 'x']
        ARRAY_DIMENSIONS = ['time', 'channel', 'z', 'y', 'x']
        chunk_tc = (
            chunk_time or data.shape[3],
            chunk_channel or data.shape[4]
        )
    elif data.ndim == 4:
        nbatch = 1
        perm = [3, 2, 1, 0]
        axes = ['t', 'z', 'y', 'x']
        ARRAY_DIMENSIONS = ['time', 'z', 'y', 'x']
        chunk_tc = (chunk_time or data.shape[3],)
    elif data.ndim == 3:
        nbatch = 0
        perm = [2, 1, 0]
        axes = ['z', 'y', 'x']
        ARRAY_DIMENSIONS = ['z', 'y', 'x']
        chunk_tc = tuple()
    elif data.ndim > 5:
        raise ValueError('Too many dimensions for conversion to nii.zarr')
    else:
        raise ValueError('Too few dimensions, this is weird')
    data = data.transpose(perm)

    # Compute image pyramid
    if label is None:
        label = jsonheader["intent"]["code"] in ("LABEL", "NEURONAMES")
    pyramid_fn = pyramid_gaussian if method[0] == 'g' else pyramid_laplacian
    data = list(_make_pyramid3d(data, nb_levels, pyramid_fn, label,
                                no_pyramid_axis))
    nb_levels = len(data)
    shapes = [d.shape[-3:] for d in data]

    # Prepare array metadata at each level
    compressor = _make_compressor(compressor, **compressor_options)
    if not isinstance(chunk, list):
        chunk = [chunk]
    chunk = [tuple(c) if isinstance(c, (list, tuple)) else (c,) for c in chunk]
    chunk = [c + c[-1:] * max(0, 3 - len(c)) + chunk_tc for c in chunk]
    chunk = [tuple(c[i] for i in perm) for c in chunk]
    chunk = chunk + chunk[-1:] * max(0, nb_levels - len(chunk))
    chunk = [{
        'chunks': c,
        'dimension_separator': r'/',
        'order': 'F',
        'dtype': jsonheader['datatype'],
        'fill_value': fill_value,
        'compressor': compressor,
    } for c in chunk]

    # Write zarr arrays
    write_multiscale(
        data, out,
        axes=axes,
        storage_options=chunk
    )

    # Write nifti header (binary)
    binheader = np.frombuffer(header.tobytes(), dtype='u1')
    out.create_dataset(
        'nifti',
        data=binheader,
        shape=[len(binheader)],
        chunks=len(binheader),
        dtype='u1',
        compressor=None,
        fill_value=None,
        dimension_separator='/',
        overwrite=True,
    )

    # Write nifti header (JSON)
    out['nifti'].attrs.update(jsonheader)

    # write xarray metadata
    for i in range(len(data)):
        out[i].attrs['_ARRAY_DIMENSIONS'] = ARRAY_DIMENSIONS

    # Ensure that OME attributes are compatible
    multiscales = out.attrs["multiscales"]
    multiscales[0]["axes"] = [
        {
            "name": "z",
            "type": "space",
            "unit": jsonheader["units"]["space"],
        },
        {
            "name": "y",
            "type": "space",
            "unit": jsonheader["units"]["space"],
        },
        {
            "name": "x",
            "type": "space",
            "unit": jsonheader["units"]["space"],
        }
    ]
    if nbatch >= 2:
        multiscales[0]["axes"].insert(0, {
            "name": "c",
            "type": "channel"
        })
    if nbatch >= 1:
        multiscales[0]["axes"].insert(0, {
            "name": "t",
            "type": "time",
            "unit": jsonheader["units"]["time"],
        })
    for n, level in enumerate(multiscales[0]["datasets"]):
        # skimage pyramid_gaussian/pyramid_laplace use
        # scipy.ndimage.zoom(..., grid_mode=True)
        # so the effective scaling is the shape ratio, and there is
        # a half voxel shift wrt to the "center of first voxel" frame
        level["coordinateTransformations"][0]["scale"] = [1.0] * nbatch + [
            (shapes[0][0]/shapes[n][0])*jsonheader["pixdim"][2],
            (shapes[0][1]/shapes[n][1])*jsonheader["pixdim"][1],
            (shapes[0][2]/shapes[n][2])*jsonheader["pixdim"][0],
        ]
        level["coordinateTransformations"].append({
            "type": "translation",
            "translation": [0.0] * nbatch + [
                (shapes[0][0]/shapes[n][0] - 1)*jsonheader["pixdim"][2]*0.5,
                (shapes[0][1]/shapes[n][1] - 1)*jsonheader["pixdim"][1]*0.5,
                (shapes[0][2]/shapes[n][2] - 1)*jsonheader["pixdim"][0]*0.5,
            ]
        })
    multiscales[0]["coordinateTransformations"] = [
        {
            "scale": [1.0] * 3,
            "type": "scale"
        }
    ]
    if nbatch >= 2:
        multiscales[0]["coordinateTransformations"][0]['scale'].insert(
            0, 1.0)
    if nbatch >= 1:
        multiscales[0]["coordinateTransformations"][0]['scale'].insert(
            0, jsonheader["pixdim"][3] or 1.0)

    out.attrs["multiscales"] = multiscales


def cli(args=None):
    """Command-line entrypoint"""
    parser = argparse.ArgumentParser(
        'nii2zarr', description='Convert nifti to nifti-zarr')
    parser.add_argument(
        'input', help='Input nifti file')
    parser.add_argument(
        'output', help='Output zarr directory')
    parser.add_argument(
        '--chunk', type=int, default=64, help='Spatial chunk size')
    parser.add_argument(
        '--unchunk-channels', action='store_true',
        help='Save all chanels in a single chunk. '
        'Unchunk if you want to display all channels as a single RGB '
        'layer in neuroglancer. Chunked by default, unless datatype is RGB.'
    )
    parser.add_argument(
        '--unchunk-time', action='store_true',
        help='Save all timepoints in a single chunk.'
        'Unchunk if you want to display all timepoints as a single RGB '
        'layer in neuroglancer. Chunked by default. '
    )
    parser.add_argument(
        '--levels', type=int, default=-1,
        help='Number of levels in the pyramid. If -1 (default), use '
             'as many levels as possible')
    parser.add_argument(
        '--method', choices=('gaussian', 'laplacian'), default='gaussian',
        help='Pyramid method')
    parser.add_argument(
        '--fill', default=None, help='Missing value')
    parser.add_argument(
        '--compressor', choices=('blosc', 'zlib'), default='blosc',
        help='Compressor')
    parser.add_argument(
        '--label', action='store_true', default=None,
        help='Segmentation volume')
    parser.add_argument(
        '--no-label', action='store_false', dest='label',
        help='Not a segmentation volume')
    parser.add_argument(
        '--no-time', action='store_true',  help='No time dimension')
    parser.add_argument(
        '--no-pyramid-axis',
        help='Thick slice axis that should not be downsampled'
    )

    args = args or sys.argv[1:]
    args = parser.parse_args(args)

    nii2zarr(
        args.input, args.output,
        chunk=args.chunk,
        chunk_channel=0 if args.unchunk_channels else 1,
        chunk_time=0 if args.unchunk_time else 1,
        nb_levels=args.levels,
        method=args.method,
        fill_value=args.fill,
        compressor=args.compressor,
        label=args.label,
        no_time=args.no_time,
        no_pyramid_axis=args.no_pyramid_axis,
    )
