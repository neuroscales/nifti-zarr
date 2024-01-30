import sys
import json
import math
import base64
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

    # Dump the binary header
    jsonheader["base64"] = base64.b64encode(header.tobytes()).decode()

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


def _make_pyramid5d(
        data5d, nb_levels, pyramid_fn=pyramid_gaussian, label=False
):
    nt, nc, *nxyz = data5d.shape
    data5d = data5d.reshape((-1, *nxyz))
    max_layer = nb_levels - 1 if nb_levels > 0 else -1

    def pyramid_values(x):
        return pyramid_fn(x, max_layer, 2, preserve_range=True)

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
                value[mask] = label[mask]
                maxprob[mask] = prob[mask]

        for level in pyramid:
            yield level

    pyramid = pyramid_labels if label else pyramid_values

    for level in zip(*map(pyramid, data5d)):
        yield np.stack(level).reshape((nt, nc) + level[0].shape)


def nii2zarr(inp, out, *,
             chunk=64,
             nb_levels=-1,
             method='gaussian',
             label=None,
             fill_value=float('nan'),
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
    chunk : [list of [list of]] int
        Chunk size, per x/y/z/t/c, per level.
        * The inner list allows different chunk sizes to be used along
          each dimension.
        * The outer list allows different chunk sizes to be used at
          different pyramid levels.
    nb_levels : int
        Number of levels in the pyramid. Default: all possible levels.
    method : {'gaussian', 'laplacian'}
        Method used to compute the pyramid.
    label : bool
        Is this is a label volume?  If `None`, guess from intent code.
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

    # Parse nifti header
    v = int(inp.header.structarr['magic'].tobytes().decode()[2])
    header = np.frombuffer(inp.header.structarr.tobytes(), count=1,
                           dtype=HEADERTYPE1 if v == 1 else HEADERTYPE2)[0]
    if header['magic'].decode() not in ('ni1', 'n+1', 'ni2', 'n+2'):
        header = header.newbyteorder()
    jsonheader = nii2json(header)

    # Fix array shape
    data = np.asarray(inp.dataobj)
    while data.ndim < 5:
        data = data[..., None]
    if data.ndim > 5:
        raise ValueError('Too many dimensions for conversion to nii.zarr')
    data = data.transpose([3, 4, 2, 1, 0])

    # Compute image pyramid
    if label is None:
        label = jsonheader["intent"]["code"] == "LABEL"
    pyramid_fn = pyramid_gaussian if method[0] == 'g' else pyramid_laplacian
    data = list(_make_pyramid5d(data, nb_levels, pyramid_fn, label))
    nb_levels = len(data)
    shapes = [d.shape[2:] for d in data]

    # Prepare array metadata at each level
    compressor = _make_compressor(compressor, **compressor_options)
    if not isinstance(chunk, (list, tuple)):
        chunk = [chunk]
    chunk = list(chunk)
    if not isinstance(chunk[0], (list, tuple)):
        chunk = [chunk]
    chunk = list(map(list, chunk))
    chunk = [c + c[-1:] * max(0, 3 - len(c)) for c in chunk]
    chunk = [c + [1] * max(0, 5 - len(c)) for c in chunk]
    chunk = [[c[i] for i in [3, 4, 2, 1, 0]] for c in chunk]
    chunk = [{
        'chunks': c,
        'dimension_separator': r'/',
        'order': 'F',
        'dtype': jsonheader['datatype'],
        'fill_value': fill_value,
        'compressor': compressor,
    } for c in chunk]
    chunk = chunk + chunk[-1:] * max(0, nb_levels - len(chunk))

    # Write zarr arrays
    write_multiscale(
        data, out,
        axes=('t', 'c', 'z', 'y', 'x'),
        storage_options=chunk
    )

    # Write nifti attributes
    out.attrs["nifti"] = jsonheader

    # write xarray metadata
    for i in range(len(data)):
        out[i].attrs['_ARRAY_DIMENSIONS'] = ['time', 'channel', 'z', 'y', 'x']

    # Ensure that OME attributes are compatible
    multiscales = out.attrs["multiscales"]
    multiscales[0]["axes"] = [
        {
            "name": "t",
            "type": "time",
            "unit": jsonheader["units"]["time"],
        },
        {
            "name": "c",
            "type": "channel"
        },
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
    for n, level in enumerate(multiscales[0]["datasets"]):
        # skimage pyramid_gaussian/pyramid_laplace use
        # scipy.ndimage.zoom(..., grid_mode=True)
        # so the effective scaling is the shape ratio, and there is
        # a half voxel shift wrt to the "center of first voxel" frame
        level["coordinateTransformations"][0]["scale"] = [
            1.0,
            1.0,
            (shapes[0][0]/shapes[n][0])*jsonheader["pixdim"][2],
            (shapes[0][1]/shapes[n][1])*jsonheader["pixdim"][1],
            (shapes[0][2]/shapes[n][2])*jsonheader["pixdim"][0],
        ]
        level["coordinateTransformations"].append({
            "type": "translation",
            "translation": [
                0.0,
                0.0,
                (shapes[0][0]/shapes[n][0] - 1)*jsonheader["pixdim"][2]*0.5,
                (shapes[0][1]/shapes[n][1] - 1)*jsonheader["pixdim"][1]*0.5,
                (shapes[0][2]/shapes[n][2] - 1)*jsonheader["pixdim"][0]*0.5,
            ]
        })
    multiscales[0]["coordinateTransformations"] = [
        {
            "scale": [
                jsonheader["pixdim"][3]
                if len(jsonheader["pixdim"]) > 3 else 1.0,
                1.0,
                1.0,
                1.0,
                1.0
            ],
            "type": "scale"
        }
    ]
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
        '--chunk', type=int, default=64, help='Chunk size')
    parser.add_argument(
        '--levels', type=int, default=-1,
        help='Number of levels in the pyramid. If -1 (default), use '
             'as many levels as possible')
    parser.add_argument(
        '--method', choices=('gaussian', 'laplacian'), default='gaussian',
        help='Pyramid method')
    parser.add_argument(
        '--fill', type=float, default=float('nan'),
        help='Missing value')
    parser.add_argument(
        '--compressor', choices=('blosc', 'zlib'), default='blosc',
        help='Compressor')

    args = args or sys.argv[1:]
    args = parser.parse_args(args)

    nii2zarr(
        args.input, args.output,
        chunk=args.chunk,
        nb_levels=args.levels,
        method=args.method,
        fill_value=args.fill,
        compressor=args.compressor,
    )
