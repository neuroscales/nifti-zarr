import NIfTI: get_sform, get_qform
import LinearAlgebra: Diagonal
import StaticArrays: SVector
import NIfTI: NIVolume


_zarr2nii_prepare_out(x) = _open_stream(x, "w")
_zarr2nii_prepare_out(x::Nothing) = x
_zarr2nii_prepare_out(x::OpenedStream) = x
_zarr2nii_prepare_inp(x) = Zarr.zopen(x)
_zarr2nii_prepare_inp(x::Zarr.ZGroup) = x


# massage inputs that have a "liberal" format
zarr2nii(inp, out=nothing; level::Integer=0) = zarr2nii(
    _zarr2nii_prepare_inp(inp),
    _zarr2nii_prepare_out(out);
    level=level
)

# version that writes out a nifti file
function zarr2nii(inp::Zarr.ZGroup, out::OpenedStream; level::Integer=0)
    vol = zarr2nii(inp; level=level)
    # We have to load it in contiguous memory, otherwise bytes are written
    # all shuffled (i.e., julia writes the contiguous bytes in the array,
    # wether or not their data are F-ordered). Since we use permutedims
    # in the main function, our data is never F-ordered so we always
    # need to make a copy before writing to disk.
    vol = NIVolume(vol.header, vol.extensions, copy(vol.raw))
    write(out.stream, vol)
    map(close, out.mine)
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
    binheader = zeros(UInt8, length(inp.arrays["nifti"]))
    copy!(binheader, inp.arrays["nifti"])
    niiheader = bytesdecode(NIfTI.NIfTI1Header, binheader)
    niiheader.magic = MAGIC1P

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
    if nbdim == 3
        perm = (3, 2, 1)
    elseif nbdim == 4
        perm = (4, 3, 2, 1)
    elseif nbdim == 5
        perm = (5, 4, 3, 1, 2)
    else
        error("Array should have 3 to 5 dimensions")
    end
    niiarray = inp.arrays[string(level)]
    niiarray = permutedims(niiarray, perm)

    # create NIVolume
    return NIVolume(niiheader, niiarray)
end
