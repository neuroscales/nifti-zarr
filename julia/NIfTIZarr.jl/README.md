# Example implementation of the nifti-zarr specification in JULIA

Convert a nifti file to a nifti-zarr storage.
```julia
import NIfTIZarr: nii2zarr
nii2zarr("path/to/nifti.nii.gz", "s3://path/to/bucket")
```

Convert a nifti-zarr storage to a nifti file.
The pyramid level can be selected with `level=L`, where 0 is the
base/finest level.
```julia
import NIfTIZarr: zarr2nii
zarr2nii("s3://path/to/bucket", "path/to/nifti.nii.gz"; level=0)
```

Load a nifti-zarr into a `NIVolume` object.
```julia
import NIfTIZarr: zarr2nii
nivol = zarr2nii("s3://path/to/bucket"; level=0)
```
