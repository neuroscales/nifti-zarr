
include("../src/NIfTIZarr.jl")
using .NIfTIZarr
import NIfTI
original_file = "/local_mount/space/megaera/1/users/kchai/nifti-zarr-julia/julia/NIfTIZarr.jl/test/data/example4d.nii.gz"
output_file = "/local_mount/space/megaera/1/users/kchai/nifti-zarr-julia/julia/NIfTIZarr.jl/test/output/example4d.nii.zarr"
# julia/NIfTIZarr.jl/src/NIfTIZarr.jl
# julia/NIfTIZarr.jl/src/NIfTIZarr.jl
# NIfTIZarr.nii2zarr("/local_mount/space/megaera/1/users/kchai/nifti-zarr-julia/julia/NIfTIZarr.jl/test/data/example4d.nii.gz","/local_mount/space/megaera/1/users/kchai/nifti-zarr-julia/julia/NIfTIZarr.jl/test/output/example4d.nii.zarr",chunk=128)
# NIfTIZarr.zarr2nii("/local_mount/space/megaera/1/users/kchai/nifti-zarr-julia/julia/NIfTIZarr.jl/test/output/example4d.nii.zarr","/local_mount/space/megaera/1/users/kchai/nifti-zarr-julia/julia/NIfTIZarr.jl/test/output/example4d.nii.gz")


NIfTIZarr.nii2zarr(original_file,output_file,chunk=32)
back =NIfTIZarr.zarr2nii(output_file)
original = NIfTI.niread(original_file)

println(original.header==back.header)
println(original.raw==back.raw)

