# Example implementation of the nifti-zarr specification in Python

Convert a nifti file to a nifti-zarr storage.

```python
from niizarr import nii2zarr
nii2zarr("path/to/nifti.nii.gz", "s3://path/to/bucket")
```

Convert a nifti-zarr storage to a nifti file. The pyramid level can be selected with `level=L`, where 0 is the
base/finest level.

```python
from niizarr import zarr2nii
zarr2nii("s3://path/to/bucket", "path/to/nifti.nii.gz", level=0)
```

Load a nifti-zarr into a `nibabel.Nifti1Image` object.

```python
from niizarr import zarr2nii
nivol = zarr2nii("s3://path/to/bucket", level=0)
```
