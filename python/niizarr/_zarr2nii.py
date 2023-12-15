import base64
import io
import zarr.hierarchy
import zarr.storage
import numpy as np
import dask.array
from nibabel import (Nifti1Image, Nifti1Header, Nifti2Image, Nifti2Header,
                     save, load)
from ._header import HEADERTYPE1, HEADERTYPE2

# If fsspec available, use fsspec
try:
    import fsspec
    open = fsspec.open
except (ImportError, ModuleNotFoundError):
    fsspec = None


def bin2nii(buffer):
    header = np.frombuffer(buffer, dtype=HEADERTYPE1, count=1)[0]
    if header['magic'].decode() not in ('ni1', 'n+1', 'nz1'):
        header = header.newbyteorder()
    if header['magic'].decode() not in ('ni1', 'n+1', 'nz1'):
        header = np.frombuffer(buffer, dtype=HEADERTYPE2, count=1)[0]
        if header['magic'].decode() not in ('ni2', 'n+2', 'nz2'):
            header = header.newbyteorder()
        if header['magic'].decode() not in ('ni2', 'n+2', 'nz2'):
            raise ValueError('Is this a nifti header?')
    return header


def zarr2nii(inp, out=None, level=0):
    """
    Convert a nifti-zarr to nifti

    Parameters
    ----------
    inp : zarr.Store or zarr.Group or path
        Output zarr object
    out : path or file_like, optional
        Path to output file. If not provided, do not write a file.
    levels : int
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

    # build structured header
    binheader = base64.b64decode(inp.attrs['nifti']['base64'])
    header = bin2nii(binheader)

    # create nibabel header (useful to convert quat 2 affine, etc)
    header = header.copy()
    magic = header['magic'].decode()
    header['magic'] = 'n+1' if magic[-1] == '1' else 'n+2'
    if magic[-1] == '1':
        NiftiHeader = Nifti1Header
        NiftiImage = Nifti1Image
    else:
        NiftiHeader = Nifti2Header
        NiftiImage = Nifti2Image
    niiheader = NiftiHeader.from_fileobj(io.BytesIO(header.tobytes()))

    # create affine at current resolution
    if level != 0:
        qform, qcode = niiheader.get_qform(coded=True)
        sform, scode = niiheader.get_sform(coded=True)
        datasets = inp.attrs['multiscales'][0]['datasets']

        xfrm0 = datasets[0]['coordinateTransformations']
        phys0 = np.eye(4)
        phys0[[0, 1, 2], [0, 1, 2]] = list(reversed(xfrm0[0]['scale'][2:]))
        if len(xfrm0) > 1:
            phys0[:3, -1] = list(reversed(xfrm0[1]['translation'][2:]))

        xfrm1 = datasets[level]['coordinateTransformations']
        phys1 = np.eye(4)
        phys1[[0, 1, 2], [0, 1, 2]] = list(reversed(xfrm1[0]['scale'][2:]))
        if len(xfrm1) > 1:
            phys1[:3, -1] = list(reversed(xfrm1[1]['translation'][2:]))

        qform = qform @ (np.linalg.inv(phys0) @ phys1)
        sform = sform @ (np.linalg.inv(phys0) @ phys1)
        niiheader.set_qform(qform, qcode)
        niiheader.set_sform(sform, scode)

    # reorder/reshape array as needed
    array = dask.array.from_zarr(inp[f'{level}'])
    array = array.transpose([4, 3, 2, 0, 1])
    for _ in range(5 - header['dim'][0].item()):
        assert array.shape[-1] == 1
        array = array[..., 0]

    # create nibabel image
    img = NiftiImage(array, None, niiheader)

    if out is not None:
        if hasattr(out, 'read'):
            img.to_stream(out)
            img = NiftiImage.from_stream(inp)
        else:
            save(img, out)
            img = load(out)

    return img
