import argparse
import io
import json
import math
import sys
import re

import numcodecs
import numpy as np
import zarr.hierarchy
import zarr.storage
from nibabel import Nifti1Image, load
from ome_zarr.writer import write_multiscale
from skimage.transform import pyramid_gaussian, pyramid_laplacian

from ._header import (
    UNITS, DTYPES, INTENTS, INTENTS_P, SLICEORDERS, XFORMS, bin2nii, get_magic_string, SYS_BYTEORDER, JNIFTI_ZARR,
    SYS_BYTEORDER_SWAPPED
)

# If fsspec available, use fsspec
try:
    import fsspec

    open = fsspec.open
except (ImportError, ModuleNotFoundError):
    fsspec = None


def nii2json(header, extensions=False):
    """
    Convert a nifti header to JSON

    Parameters
    ----------
    header : np.array[HEADERTYPE1 or HEADERTYPE2]
        Nifti header in binary form
    extensions: bool, optional
        If extensions present in this nifti file
    Returns
    -------
    header : dict
        Nifti header in JSON form

    """
    header = header.copy()

    ndim = header["dim"][0].item()
    intent_code = INTENTS[header["intent_code"].item()]

    intent_param = header["intent_p"][:INTENTS_P[intent_code]].tolist()
    quatern = header["quatern"].tolist()
    qoffset = header["qoffset"].tolist()

    jsonheader = {
        # Strip control characters
        "NIIFormat": get_magic_string(header),
        "Dim": header["dim"][1:1 + ndim].tolist(),
        "VoxelSize": header["pixdim"][1:1 + ndim].tolist(),
        "Unit": {
            "L": UNITS[(header["xyzt_units"] & 0x07).item()],
            "T": UNITS[(header["xyzt_units"] & 0x38).item()],
        },
        "DataType": DTYPES[header["datatype"].item()],
        "DimInfo": {
            "Freq": (header["dim_info"] & 0x03).item(),
            "Phase": ((header["dim_info"] >> 2) & 0x03).item(),
            "Slice": ((header["dim_info"] >> 4) & 0x03).item(),
        },
        "Intent": intent_code,
        "Param1": intent_param[0] if len(intent_param) > 0 else None,
        "Param2": intent_param[1] if len(intent_param) > 1 else None,
        "Param3": intent_param[2] if len(intent_param) > 2 else None,
        "Name": header["intent_name"].tobytes().decode(),
        "ScaleSlope": header["scl_slope"].item(),
        "ScaleOffset": header["scl_inter"].item(),
        "FirstSliceID": header["slice_start"].item(),
        "LastSliceID": header["slice_end"].item(),
        "SliceType": SLICEORDERS[header["slice_code"].item()],
        "SliceTime": header["slice_duration"].item(),
        "MinIntensity": header["cal_min"].item(),
        "MaxIntensity": header["cal_max"].item(),
        "TimeOffset": header["toffset"].item(),
        "Description": header["descrip"].tobytes().decode(),
        "AuxFile": header["aux_file"].tobytes().decode(),
        "QForm": XFORMS[header["qform_code"].item()],
        "Quatern": {
            "b": quatern[0],
            "c": quatern[1],
            "d": quatern[2],
        },
        "QuaternOffset": {
            "x": qoffset[0],
            "y": qoffset[1],
            "z": qoffset[2],
        },
        # TODO: change this after JNIfTI changed
        "Orientation": {
            "x": "r" if header["pixdim"][0].item() == 0 else "l",
            "y": "a",
            "z": "s",
        },
        "SForm": XFORMS[header["sform_code"].item()],
        "Affine": header["sform"].tolist(),
        "NIFTIExtension": [1 if extensions else 0] + [0, 0, 0],
    }
    if not math.isfinite(jsonheader["ScaleSlope"]):
        jsonheader["ScaleSlope"] = 0.0
    if not math.isfinite(jsonheader["ScaleOffset"]):
        jsonheader["ScaleOffset"] = 0.0

    # Remove control characters
    for k,v in jsonheader.items():
        if isinstance(v, str):
            jsonheader[k] = re.sub(r'[\n\r\t\00]*', '', v)

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
        'x': 2,
        'y': 1,
        'z': 0,
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
    no_pyramid_axis : {'x', 'y', 'z'}
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

    header, byteorder_swapped = bin2nii(inp.header.structarr.tobytes(), True)

    jsonheader = nii2json(header, len(inp.header.extensions) != 0)

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
        label = jsonheader['Intent'] in ("label", "neuronames")
    pyramid_fn = pyramid_gaussian if method[0] == 'g' else pyramid_laplacian
    data = list(_make_pyramid3d(data, nb_levels, pyramid_fn, label,
                                no_pyramid_axis))
    nb_levels = len(data)
    shapes = [d.shape[-3:] for d in data]

    # Fix data type
    # If nifti was swapped when loading it, we want to swapped it back to make it as same as before
    byteorder = SYS_BYTEORDER_SWAPPED if byteorder_swapped else SYS_BYTEORDER
    data_type = JNIFTI_ZARR[jsonheader['DataType']]
    if isinstance(data_type, tuple):
        data_type = [
            [field, '|' + dtype] for field, dtype in data_type
        ]
    elif data_type.endswith('1'):
        data_type = '|' + data_type
    else:
        data_type = byteorder + data_type

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
        'dtype': data_type,
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
    header_data = header.tobytes().newbyteorder() if byteorder_swapped else header.tobytes()
    bin_data = [np.frombuffer(header_data, dtype='u1')]

    if inp.header.extensions:
        extension_stream = io.BytesIO()
        inp.header.extensions.write_to(extension_stream,
                                       byteswap=byteorder_swapped)
        bin_data.append(np.frombuffer(extension_stream.getvalue(), dtype=np.uint8))

    # Concatenate the final binary data
    bin_data = np.concatenate(bin_data)

    out.create_dataset(
        'nifti',
        data=bin_data,
        shape=[len(bin_data)],
        chunks=len(bin_data),
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
            "unit": JNIFTI_ZARR[jsonheader["Unit"]["L"]],
        },
        {
            "name": "y",
            "type": "space",
            "unit": JNIFTI_ZARR[jsonheader["Unit"]["L"]],
        },
        {
            "name": "x",
            "type": "space",
            "unit": JNIFTI_ZARR[jsonheader["Unit"]["L"]],
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
            "unit": JNIFTI_ZARR[jsonheader["Unit"]["T"]],
        })
    for n, level in enumerate(multiscales[0]["datasets"]):
        # skimage pyramid_gaussian/pyramid_laplace use
        # scipy.ndimage.zoom(..., grid_mode=True)
        # so the effective scaling is the shape ratio, and there is
        # a half voxel shift wrt to the "center of first voxel" frame
        level["coordinateTransformations"][0]["scale"] = [1.0] * nbatch + [
            (shapes[0][0] / shapes[n][0]) * jsonheader["VoxelSize"][2],
            (shapes[0][1] / shapes[n][1]) * jsonheader["VoxelSize"][1],
            (shapes[0][2] / shapes[n][2]) * jsonheader["VoxelSize"][0],
        ]
        level["coordinateTransformations"].append({
            "type": "translation",
            "translation": [0.0] * nbatch + [
                (shapes[0][0] / shapes[n][0] - 1) * jsonheader["VoxelSize"][2] * 0.5,
                (shapes[0][1] / shapes[n][1] - 1) * jsonheader["VoxelSize"][1] * 0.5,
                (shapes[0][2] / shapes[n][2] - 1) * jsonheader["VoxelSize"][0] * 0.5,
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
            0, jsonheader["VoxelSize"][3] or 1.0)

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
        '--no-time', action='store_true', help='No time dimension')
    parser.add_argument(
        '--no-pyramid-axis', choices=('x', 'y', 'z'),
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
