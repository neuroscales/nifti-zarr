import argparse
import io
import json
import math
import re
import sys

import numcodecs
import numpy as np
import zarr
import zarr.storage
import zarr.codecs
from nibabel import Nifti1Image, load
from skimage.transform import pyramid_gaussian, pyramid_laplacian
from ._header import (
    UNITS, DTYPES, INTENTS, INTENTS_P, SLICEORDERS, XFORMS,
    bin2nii, get_magic_string, SYS_BYTEORDER, JNIFTI_ZARR,
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
    nii_version = 1 if header["sizeof_hdr"].item() == 348 else 2
    jsonheader = {
        "NIIHeaderSize": header["sizeof_hdr"].item(),
        "DimInfo": {
            "Freq": (header["dim_info"] & 0x03).item(),
            "Phase": ((header["dim_info"] >> 2) & 0x03).item(),
            "Slice": ((header["dim_info"] >> 4) & 0x03).item(),
        },
        "Dim": header["dim"][1:1 + ndim].tolist(),
        "Param1": intent_param[0] if len(intent_param) > 0 else None,
        "Param2": intent_param[1] if len(intent_param) > 1 else None,
        "Param3": intent_param[2] if len(intent_param) > 2 else None,
        "Intent": intent_code,
        "DataType": DTYPES[header["datatype"].item()],
        "BitDepth": header["bitpix"].item(),
        "FirstSliceID": header["slice_start"].item(),
        "VoxelSize": header["pixdim"][1:1 + ndim].tolist(),
        "Orientation": {
            "x": "r" if header["pixdim"][0].item() == 0 else "l",
            "y": "a",
            "z": "s",
        },
        "NIIByteOffset": header["vox_offset"].item(),
        "ScaleSlope": header["scl_slope"].item(),
        "ScaleOffset": header["scl_inter"].item(),
        "LastSliceID": header["slice_end"].item(),
        "SliceType": SLICEORDERS[header["slice_code"].item()],
        "Unit": {
            "L": UNITS[(header["xyzt_units"] & 0x07).item()],
            "T": UNITS[(header["xyzt_units"] & 0x38).item()],
        },
        "MaxIntensity": header["cal_max"].item(),
        "MinIntensity": header["cal_min"].item(),
        "SliceTime": header["slice_duration"].item(),
        "TimeOffset": header["toffset"].item(),
        "Description": header["descrip"].tobytes().decode(),
        "AuxFile": header["aux_file"].tobytes().decode(),
        "QForm": XFORMS[header["qform_code"].item()],
        "SForm": XFORMS[header["sform_code"].item()],
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
        "Affine": header["sform"].tolist(),
        "Name": header["intent_name"].tobytes().decode(),
        # Strip control characters
        "NIIFormat": get_magic_string(header),
        "NIFTIExtension": [1 if extensions else 0] + [0, 0, 0],
    }
    if not math.isfinite(jsonheader["ScaleSlope"]):
        jsonheader["ScaleSlope"] = 0.0
    if not math.isfinite(jsonheader["ScaleOffset"]):
        jsonheader["ScaleOffset"] = 0.0

    if nii_version == 1:
        unused_fields = {
            "A75DataTypeName": header["datatype"].tobytes().decode(),
            "A75DBName": header["db_name"].tobytes().decode(),
            "A75Extends": header["extents"].item(),
            "A75SessionError": header["session_error"].item(),
            "A75Regular": header["regular"].item(),
            "A75GlobalMax": header["glmax"].item(),
            "A75GlobalMin": header["glmin"].item(),
        }
    else:
        unused_fields = {
            "A75DataTypeName": "",
            "A75DBName": "",
            "A75Extends": 0,
            "A75SessionError": 0,
            "A75Regular": 0,
            "A75GlobalMax": 0,
            "A75GlobalMin": 0,
        }
    jsonheader.update(unused_fields)
    # Remove control characters
    for k, v in jsonheader.items():
        if isinstance(v, str):
            jsonheader[k] = re.sub(r'[\n\r\t\00]*', '', v)

    # Check that the dictionary is serializable
    json.dumps(jsonheader)

    return jsonheader


def _make_compressor(name, zarr_version, **prm):
    if not isinstance(name, str):
        return name
    name = name.lower()
    if name == 'blosc':
        if zarr_version == 3:
            Compressor = zarr.codecs.BloscCodec
        elif zarr_version == 2:
            Compressor = numcodecs.Blosc
    elif name == 'zlib':
        if zarr_version == 3:
            Compressor = zarr.codecs.ZstdCodec
        elif zarr_version == 2:
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
    max_layer = nb_levels - 1

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


def write_ome_metadata(
        omz: zarr.Group,
        axes: list[str],
        space_scale: float | list[float] = 1,
        time_scale: float = 1,
        space_unit: str = "micrometer",
        time_unit: str = "second",
        name: str = "",
        pyramid_aligns: str | int | list[str | int] = 2,
        levels: int | None = None,
        no_pool: int | None = None,
        multiscales_type: str = "",
        ome_version: str = "0.4"
) -> None:
    """
    Write OME metadata into Zarr.

    Parameters
    ----------
    path : str | PathLike
        Path to parent Zarr.
    axes : list[str]
        Name of each dimension, in Zarr order (t, c, z, y, x)
    space_scale : float | list[float]
        Finest-level voxel size, in Zarr order (z, y, x)
    time_scale : float
        Time scale
    space_unit : str
        Unit of spatial scale (assumed identical across dimensions)
    space_time : str
        Unit of time scale
    name : str
        Name attribute
    pyramid_aligns : float | list[float] | {"center", "edge"}
        Whether the pyramid construction aligns the edges or the centers
        of the corner voxels. If a (list of) number, assume that a moving
        window of that size was used.
    levels : int
        Number of existing levels. Default: find out automatically.
    zarr_version : {2, 3} | None
        Zarr version. If `None`, guess from existing zarr array.

    """
    # Read shape at each pyramid level
    shapes = []
    level = 0
    while True:
        if levels is not None and level > levels:
            break

        if str(level) not in omz.keys():
            levels = level
            break
        shapes += [
            omz[str(level)].shape,
        ]
        level += 1

    axis_to_type = {
        "x": "space",
        "y": "space",
        "z": "space",
        "t": "time",
        "c": "channel",
    }

    # Number of spatial (s), batch (b) and total (n) dimensions
    ndim = len(axes)
    sdim = sum(axis_to_type[axis] == "space" for axis in axes)
    bdim = ndim - sdim

    if isinstance(pyramid_aligns, (int, str)):
        pyramid_aligns = [pyramid_aligns]
    pyramid_aligns = list(pyramid_aligns)
    if len(pyramid_aligns) < sdim:
        repeat = pyramid_aligns[:1] * (sdim - len(pyramid_aligns))
        pyramid_aligns = repeat + pyramid_aligns
    pyramid_aligns = pyramid_aligns[-sdim:]

    if isinstance(space_scale, (int, float)):
        space_scale = [space_scale]
    space_scale = list(space_scale)
    if len(space_scale) < sdim:
        repeat = space_scale[:1] * (sdim - len(space_scale))
        space_scale = repeat + space_scale
    space_scale = space_scale[-sdim:]

    multiscales = [
        {
            "version": ome_version,
            "axes": [
                {
                    "name": axis,
                    "type": axis_to_type[axis],
                }
                if axis_to_type[axis] == "channel"
                else {
                    "name": axis,
                    "type": axis_to_type[axis],
                    "unit": (
                        space_unit
                        if axis_to_type[axis] == "space"
                        else time_unit
                        if axis_to_type[axis] == "time"
                        else None
                    ),
                }
                for axis in axes
            ],
            "datasets": [],
            "type": "median window " + "x".join(["2"] * sdim)
            if not multiscales_type
            else multiscales_type,
            "name": name,
        }
    ]

    shape0 = shapes[0]
    for n in range(len(shapes)):
        shape = shapes[n]
        multiscales[0]["datasets"].append({})
        level = multiscales[0]["datasets"][-1]
        level["path"] = str(n)

        scale = [1] * bdim + [
            (
                pyramid_aligns[i] ** n
                if not isinstance(pyramid_aligns[i], str)
                else (shape0[bdim + i] / shape[bdim + i])
                if pyramid_aligns[i][0].lower() == "e"
                else ((shape0[bdim + i] - 1) / (shape[bdim + i] - 1))
            )
            * space_scale[i]
            if i != no_pool
            else space_scale[i]
            for i in range(sdim)
        ]
        translation = [0] * bdim + [
            (
                pyramid_aligns[i] ** n - 1
                if not isinstance(pyramid_aligns[i], str)
                else (shape0[bdim + i] / shape[bdim + i]) - 1
                if pyramid_aligns[i][0].lower() == "e"
                else 0
            )
            * 0.5
            * space_scale[i]
            if i != no_pool
            else 0
            for i in range(sdim)
        ]

        level["coordinateTransformations"] = [
            {
                "type": "scale",
                "scale": scale,
            },
            {
                "type": "translation",
                "translation": translation,
            },
        ]

    scale = [1.0] * ndim
    if "t" in axes:
        scale[axes.index("t")] = time_scale
    multiscales[0]["coordinateTransformations"] = [{"scale": scale, "type": "scale"}]

    multiscales[0]["version"] = ome_version
    if ome_version == "0.4":
        omz.attrs["multiscales"] = multiscales
    elif ome_version == "0.5":
        omz.attrs["multiscales"] = multiscales
        omz.attrs["ome"] = {"multiscales": multiscales}
    else:
        raise Exception("Unsupported ome version")


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
             compressor_options={},
             zarr_version=2,
             ome_version="0.4",
             ):
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
        Number of levels in the pyramid.
        If -1, make all possible levels until the level can be fit into
        one chunk.
        Default: -1
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
    zarr_version : {2, 3}
        Zarr format version
    ome_version : {0.4, 0.5}
        OME-Zarr version
    """
    # Open nifti image with nibabel
    if not isinstance(inp, Nifti1Image):
        if hasattr(inp, 'read'):
            inp = Nifti1Image.from_stream(inp)
        else:
            inp = load(inp)

    # Open zarr group
    if not isinstance(out, zarr.Group):
        if not isinstance(out, zarr.storage.StoreLike):
            if fsspec:
                out = zarr.storage.FsspecStore(out)
            else:
                out = zarr.storage.LocalStore(out)
        out = zarr.group(store=out, overwrite=True, zarr_version=zarr_version)

    if no_time and len(inp.shape) > 3:
        inp = Nifti1Image(inp.dataobj[:, :, :, None], inp.affine, inp.header)
    # nibabel consumde these two values
    if hasattr(inp.dataobj, "_slope") and hasattr(inp.dataobj, "_inter"):
        inp.header.set_slope_inter(inp.dataobj._slope, inp.dataobj._inter)

    header = bin2nii(inp.header.structarr.tobytes())
    byteorder_swapped = inp.header.endianness != SYS_BYTEORDER
    jsonheader = nii2json(header, len(inp.header.extensions) != 0)

    if hasattr(inp.dataobj, "get_unscaled"):
        data = np.asarray(inp.dataobj.get_unscaled())
    else:
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

    # if chunk is a list, we use the first tuple,
    # otherwise the logic might be too complicated
    chunksize = np.array((chunk,) * 3 if isinstance(chunk, int) else chunk[0])
    nxyz = np.array(data.shape[-3:])

    if nb_levels == -1:
        default_nb_levels = int(np.ceil(np.log2(np.max(nxyz / chunksize)))) + 1
        default_nb_levels = max(default_nb_levels, 1)
        nb_levels = default_nb_levels

    data = list(_make_pyramid3d(data, nb_levels, pyramid_fn, label,
                                no_pyramid_axis))

    shapes = [d.shape[-3:] for d in data]

    # Fix data type
    # If nifti was swapped when loading it, we want to swapped it back
    # to make it as same as before
    byteorder = SYS_BYTEORDER_SWAPPED if byteorder_swapped else SYS_BYTEORDER
    data_type = JNIFTI_ZARR[jsonheader['DataType']]
    if isinstance(data_type, tuple):
        data_type = [
            (field, '|' + dtype) for field, dtype in data_type
        ]
    elif data_type.endswith('1'):
        data_type = '|' + data_type
    else:
        data_type = byteorder + data_type

    # Prepare array metadata at each level
    compressor = _make_compressor(compressor, zarr_version=zarr_version,
                                  **compressor_options)
    if not isinstance(chunk, list):
        chunk = [chunk]
    chunk = [tuple(c) if isinstance(c, (list, tuple)) else (c,) for c in chunk]
    chunk = [c + c[-1:] * max(0, 3 - len(c)) + chunk_tc for c in chunk]
    chunk = [tuple(c[i] for i in perm) for c in chunk]
    chunk = chunk + chunk[-1:] * max(0, nb_levels - len(chunk))
    chunk = [{
        'chunks': c,
        # 'dimension_separator': r'/',
        'order': 'C',
        'dtype': data_type,
        'fill_value': fill_value,
        'compressors': compressor,
    } for c in chunk]

    # Write zarr arrays
    # write_multiscale(
    #     data, out,
    #     axes=axes,
    #     storage_options=chunk
    # )
    for i, d in enumerate(data):
        out.create_array(str(i), shape=d.shape, **chunk[i])
        out[str(i)][:] = d
    # Write nifti header (binary)
    stream = io.BytesIO()
    inp.header.write_to(stream)
    bin_data = np.frombuffer(stream.getvalue(), dtype=np.uint8)
    # Remove the extension flag if there are no extensions
    if len(inp.header.extensions) == 0:
        bin_data = bin_data[:-4]

    out.create_array(
        'nifti',
        shape=[len(bin_data)],
        chunks=len(bin_data),
        dtype='u1',
        compressors=None,
        fill_value=None,
        # dimension_separator='/',
        overwrite=True,
    )
    out['nifti'][:] = bin_data
    # Write nifti header (JSON)
    out['nifti'].attrs.update(jsonheader)

    # write xarray metadata
    for i in range(len(data)):
        out[str(i)].attrs['_ARRAY_DIMENSIONS'] = ARRAY_DIMENSIONS

    write_ome_metadata(out,
                       axes=axes,
                       space_scale=[jsonheader["VoxelSize"][2],
                                    jsonheader["VoxelSize"][1],
                                    jsonheader["VoxelSize"][0]],
                       time_scale=jsonheader["VoxelSize"][3] if nbatch >= 1 else 1.0,
                       space_unit=JNIFTI_ZARR[jsonheader["Unit"]["L"]],
                       time_unit=JNIFTI_ZARR[jsonheader["Unit"]["T"]],
                       ome_version=ome_version
                       )
    return


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
             'layer in neuroglancer. '
             'Chunked by default, unless datatype is RGB.'
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
    parser.add_argument(
        '--zarr-version', type=int, default=2
    )
    parser.add_argument(
        '--ome-version', type=str, default="0.4"
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
        zarr_version=args.zarr_version,
        ome_version=args.ome_version,
    )
