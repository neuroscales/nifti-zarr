import argparse
import io
import sys
import warnings

import dask.array
import numpy as np
import zarr.hierarchy
import zarr.storage
from nibabel import (Nifti1Image, Nifti1Header, Nifti2Image, Nifti2Header,
                     save, load)
from nibabel.nifti1 import Nifti1Extension, Nifti1Extensions
from ._header import bin2nii, NIFTI_1_HEADER_SIZE, NIFTI_2_HEADER_SIZE

# If fsspec available, use fsspec
try:
    import fsspec
    open = fsspec.open
except (ImportError, ModuleNotFoundError):
    fsspec = None

def extract_extension(chunk, index=0):
    sections = []
    chunk_len = len(chunk)

    while index < chunk_len:
        # TODO: use numpy to load two int and fix endian?
        size = int.from_bytes(chunk[index:index + 4], byteorder='big')
        code = int.from_bytes(chunk[index + 4:index + 8], byteorder='big')

        content = chunk[index + 8:index + size]
        # strip redundant \0. See: https://github.com/nipy/nibabel/blob/83eaf0b55be9e9079bf9ad64975b71c22523f5f0/nibabel/nifti1.py#L630
        content = content.rstrip(b'\x00')
        sections.append(Nifti1Extension(code, content))

        index += size

    return sections


def zarr2nii(inp, out=None, level=0):
    """
    Convert a nifti-zarr to nifti

    Parameters
    ----------
    inp : zarr.Store or zarr.Group or path
        Output zarr object
    out : path or file_like, optional
        Path to output file. If not provided, do not write a file.
    level : int
        Pyramid level to extract

    Returns
    -------
    out : nib.Nifti1Image
        Mapped output file _or_ Nifti object whose dataobj is a dask array
    """

    if not isinstance(inp, zarr.hierarchy.Group):
        if not isinstance(inp, zarr.storage.Store):
            if fsspec:
                inp = zarr.storage.FSStore(inp)
            else:
                inp = zarr.storage.DirectoryStore(inp)
        inp = zarr.group(store=inp)

    # read binary header
    header = bin2nii(np.asarray(inp['nifti']).tobytes())

    # create nibabel header (useful to convert quat 2 affine, etc)
    if header['sizeof_hdr'] == NIFTI_1_HEADER_SIZE:
        NiftiHeader = Nifti1Header
        NiftiImage = Nifti1Image
    elif header['sizeof_hdr'] == NIFTI_2_HEADER_SIZE:
        NiftiHeader = Nifti2Header
        NiftiImage = Nifti2Image
    else:
        raise ValueError(f"sizeof_hdr {header['sizeof_hdr']} does not match any Nifti header specification")

    niiheader = NiftiHeader.from_fileobj(io.BytesIO(header.tobytes()),
                                         check=False)

    # create affine at current resolution
    if level != 0:
        qform, qcode = niiheader.get_qform(coded=True)
        sform, scode = niiheader.get_sform(coded=True)
        datasets = inp.attrs['multiscales'][0]['datasets']

        xfrm0 = datasets[0]['coordinateTransformations']
        phys0 = np.eye(4)
        phys0[[0, 1, 2], [0, 1, 2]] = list(reversed(xfrm0[0]['scale'][-3:]))
        if len(xfrm0) > 1:
            phys0[:3, -1] = list(reversed(xfrm0[1]['translation'][-3:]))

        xfrm1 = datasets[level]['coordinateTransformations']
        phys1 = np.eye(4)
        phys1[[0, 1, 2], [0, 1, 2]] = list(reversed(xfrm1[0]['scale'][-3:]))
        if len(xfrm1) > 1:
            phys1[:3, -1] = list(reversed(xfrm1[1]['translation'][-3:]))

        qform = qform @ (np.linalg.inv(phys0) @ phys1)
        sform = sform @ (np.linalg.inv(phys0) @ phys1)
        niiheader.set_qform(qform, qcode)
        niiheader.set_sform(sform, scode)

    # reorder/reshape array as needed
    array = dask.array.from_zarr(inp[f'{level}'])

    if array.ndim == 5:
        array = array.transpose([4, 3, 2, 0, 1])
        assert inp[level].attrs['_ARRAY_DIMENSIONS'] == ['time', 'channel', 'z', 'y', 'x']
    elif array.ndim == 4:
        array = array.transpose([3, 2, 1, 0])
        assert inp[level].attrs['_ARRAY_DIMENSIONS'] == ['time', 'z', 'y', 'x']
    elif array.ndim == 3:
        array = array.transpose([2, 1, 0])
        assert inp[level].attrs['_ARRAY_DIMENSIONS'] == ['z', 'y', 'x']

    # create nibabel image
    img = NiftiImage(array, None, niiheader)

    extension_size = len(inp['nifti'])-header['sizeof_hdr']
    if extension_size > 0:
        try:
            file_obj = io.BytesIO(np.asarray(inp['nifti']).tobytes()[header['sizeof_hdr']:])
            img.header.extensions = Nifti1Extensions.from_fileobj(file_obj,extension_size,False)
            # extensions = extract_extension(np.asarray(inp['nifti']).tobytes(), header['sizeof_hdr'])
            # img.header.extensions += extensions
        except Exception:
            warnings.warn("Failed to load extensions")

    if out is not None:
        if hasattr(out, 'read'):
            img.to_stream(out)
            img = NiftiImage.from_stream(inp)
        else:
            save(img, out)
            img = load(out)



    return img


def cli(args=None):
    """Command-line entrypoint"""
    parser = argparse.ArgumentParser(
        'zarr2nii', description='Convert nifti to nifti-zarr')
    parser.add_argument(
        'input', help='Input zarr directory')
    parser.add_argument(
        'output', help='Output nifti file')
    parser.add_argument(
        '--level', type=int, default=0,
        help='Pyramid level to extract (default: 0 = coarsest)')

    args = args or sys.argv[1:]
    args = parser.parse_args(args)
    zarr2nii(args.input, args.output, args.level)
