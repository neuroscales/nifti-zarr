import NIfTI: get_sform, get_qform
import LinearAlgebra: Diagonal
import StaticArrays: SVector
import NIfTI: NIVolume


_zarr2nii_prepare_out(x::AbstractString) = _stream(x, "w")
_zarr2nii_prepare_out(x::IO) = x
_zarr2nii_prepare_out(x::Nothing) = x
_zarr2nii_prepare_inp(x) = Zarr.zopen(x)
_zarr2nii_prepare_inp(x::Zarr.ZGroup) = x


# massage inputs that have a "liberal" format
zarr2nii(inp, out=nothing; level::Integer=0) = zarr2nii(
    _zarr2nii_prepare_inp(inp),
    _zarr2nii_prepare_out(out);
    level=level
)

# version that writes out a nifti file
function zarr2nii(inp::Zarr.ZGroup, out::IO; level::Integer=0)
    vol = zarr2nii(inp; level=level)
    write(out, vol)
    return vol
end

"""
    zarr2nii(inp[, out]; level=0)

Convert a nifti-zarr to nifti.

This function effectively wraps a nifti-zarr into a `NIfTI.NIVolume` object.

- The input is a Zarr store or group (or a path).
- The output is a writable IO (or a path).
  If not provided, the output nifti object is not written to disk.
- The pyramid level to extract can be provided. Level 1 is the finest level.

This function always returns a `NIfTI.NIVolume`.
"""
function zarr2nii(inp::Zarr.ZGroup, ::Nothing=nothing; level::Integer=0)

    # build structured header
    header = base64decode(inp.attrs["nifti"]["base64"], Nifti1Header)

    # create nibabel header (useful to convert quat 2 affine, etc)
    niiheader = convert(NIfTI.NIfTI1Header, header)

    # create affine at current resolution
    if level != 0
        qform = get_qform(niiheader)
        sform = get_sform(niiheader)
        datasets = inp.attrs["multiscales"][1]["datasets"]

        xfrm0 = datasets[1]["coordinateTransformations"]
        diag0 = reverse(xfrm0[1]["scale"][3:end])
        diag0 = Base.convert(Vector{Float64}, diag0)
        if length(xfrm0) > 1
            shift0 = reverse(xfrm0[2]["translation"][3:end])
        else
            shift0 = [0, 0, 0]
        end
        shift0 = Base.convert(Vector{Float64}, shift0)
        phys0 = hcat(Diagonal(diag0), shift0)
        phys0 = vcat(phys0, [0, 0, 0, 1]')

        xfrm1 = datasets[level+1]["coordinateTransformations"]
        diag1 = reverse(xfrm1[1]["scale"][3:end])
        diag1 = Base.convert(Vector{Float64}, diag1)
        if length(xfrm1) > 1
            shift1 = reverse(xfrm1[2]["translation"][3:end])
        else
            shift1 = [0, 0, 0]
        end
        shift1 = Base.convert(Vector{Float64}, shift1)
        phys1 = hcat(Diagonal(diag1), shift1)
        phys1 = vcat(phys1, [0, 0, 0, 1]')

        qform = qform * (inv(phys0) * phys1)
        sform = sform * (inv(phys0) * phys1)
        set_qform!(niiheader, qform)
        set_sform!(niiheader, sform)
    end

    # reorder/reshape array as needed
    nbdim = niiheader.dim[1]
    niiarray = inp.arrays[string(level)]
    niiarray = permutedims(niiarray, (5, 4, 3, 1, 2))
    slicer = (Tuple(Colon() for _ in 1:nbdim)...,
              Tuple(1 for _ in 1:(5-nbdim))...)
    niiarray = niiarray[slicer...]

    # create NIVolume
    return NIVolume(niiheader, niiarray)
end
