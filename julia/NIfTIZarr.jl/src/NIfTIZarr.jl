module NIfTIZarr

include("header.jl")
include("affines.jl")
include("nii2zarr.jl")
include("zarr2nii.jl")

export nii2json, nii2zarr

end  # module
