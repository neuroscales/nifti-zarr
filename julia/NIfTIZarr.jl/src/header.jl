import Base: @kwdef, read, write
import EnumX: @enumx
import Formatting: printfmt
import Base64: base64encode, Base64DecodePipe
import NIfTI
# import GZip
import CodecZlib: GzipCompressorStream, GzipDecompressorStream


"""A generic IO type"""
IOType = Union{IO,AbstractString}

struct OpenedStream{T<:IO}
    stream::T
    mine::Vector{<:IO}
end

abstract type NiftiHeader end


function Base.show(x::NiftiHeader)
    fieldlength = maximum( map(f->length(string(f)), fieldnames(typeof(x))) ) + 1
    print("$(typeof(x))()\n")
    for field in fieldnames(typeof(x))
        printfmt("    {:$(fieldlength)s} => $(reprfield(getproperty(x, field)))\n", field)
    end
end

"""
    ztuple(n::Integer)

Create a NTuple of zeros
"""
ztuple(n) = ntuple(x->0, n)

# Valid magic strings
MAGIC1 = UInt8.(('n', 'i', '1', '\0'))   # two files
MAGIC1P = UInt8.(('n', '+', '1', '\0'))  # single file
MAGIC1Z = UInt8.(('n', 'z', '1', '\0'))  # zarr folder
MAGIC2 = UInt8.(('n', 'i', '1', '\0', '\0', '\0', '\0', '\0'))   # two files
MAGIC2P = UInt8.(('n', '+', '1', '\0', '\0', '\0', '\0', '\0'))  # single file
MAGIC2Z = UInt8.(('n', 'z', '1', '\0', '\0', '\0', '\0', '\0'))  # zarr folder

"""
A NIfTI1 header in raw binary form

https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
"""
@kwdef mutable struct Nifti1Header <: NiftiHeader
    sizeof_hdr::Int32 = 348                     # MUST be 348
    data_type::NTuple{10,UInt8} = ztuple(10)    # ++UNUSED++
    db_name::NTuple{18,UInt8} = ztuple(18)      # ++UNUSED++
    extents::Int32 = 16384                      # ++UNUSED++
    session_error::Int16 = 0                    # ++UNUSED++
    regular::UInt8 = UInt8('r')                 # ++UNUSED++
    dim_info::UInt8 = 0                         # MRI slice ordering.
    dim::NTuple{8,Int16} = ztuple(8)            # Data array dimensions.
    intent_p::NTuple{3,Float32} = ztuple(3)     # Intent parameters.
    intent_code::Int16 = 0                      # NIFTI_INTENT_* code.
    datatype::Int16 = 0                         # Defines data type!
    bitpix::Int16 = 0                           # Number bits/voxel.
    slice_start::Int16 = 0                      # First slice index.
    pixdim::NTuple{8,Float32} = ztuple(8)       # Grid spacings.
    vox_offset::Float32 = 0                     # Offset into .nii file
    scl_slope::Float32 = 0                      # Data scaling: slope.
    scl_inter::Float32 = 0                      # Data scaling: offset.
    slice_end::Int16 = 0                        # Last slice index.
    slice_code::UInt8 = 0                       # Slice timing order.
    xyzt_units::UInt8 = 0                       # Units of pixdim[1..4]
    cal_max::Float32 = 0                        # Max display intensity
    cal_min::Float32 = 0                        # Min display intensity
    slice_duration::Float32 = 0                 # Time for 1 slice.
    toffset::Float32 = 0                        # Time axis shift.
    glmax::Int32 = 0                            # ++UNUSED++
    glmin::Int32 = 0                            # ++UNUSED++
    descrip::NTuple{80,UInt8} = ztuple(80)      # any text you like.
    aux_file::NTuple{24,UInt8} = ztuple(24)     # auxiliary filename.
    qform_code::Int16 = 0                       # NIFTI_XFORM_* code.
    sform_code::Int16 = 0                       # NIFTI_XFORM_* code.
    quatern::NTuple{3,Float32} = ztuple(3)      # Quaternion b/c/d param.
    qoffset::NTuple{3,Float32} = ztuple(3)      # Quaternion x/y/z param.
    srow_x::NTuple{4,Float32} = ztuple(4)       # 1st row affine transform.
    srow_y::NTuple{4,Float32} = ztuple(4)       # 2nd row affine transform.
    srow_z::NTuple{4,Float32} = ztuple(4)       # 3rd row affine transform.
    intent_name::NTuple{16,UInt8} = ztuple(16)  # 'name' or meaning of data.
    magic::NTuple{4,UInt8} = MAGIC1             # MUST be "ni1\0" or "n+1\0".
end                                             # 348 bytes total

"""
A NIfTI2 header in raw binary form

https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h
"""
@kwdef mutable struct Nifti2Header <: NiftiHeader
    sizeof_hdr::Int32 = 540                     # MUST be 540
    magic::NTuple{8,UInt8} = MAGIC2             # MUST be valid signature.
    datatype::Int16 = 0                         # Defines data type!
    bitpix::Int16 = 0                           # Number bits/voxel.
    dim::NTuple{8,Int64} = ztuple(8)            # Data array dimensions.
    intent_p::NTuple{3,Float64} = ztuple(3)     # Intent parameters.
    pixdim::NTuple{8,Float64} = ztuple(8)       # Grid spacings.
    vox_offset::Int64 = 0                       # Offset into .nii file
    scl_slope::Float64 = 0                      # Data scaling: slope.
    scl_inter::Float64 = 0                      # Data scaling: offset.
    cal_max::Float64 = 0                        # Max display intensity
    cal_min::Float64 = 0                        # Min display intensity
    slice_duration::Float64 = 0                 # Time for 1 slice.
    toffset::Float64 = 0                        # Time axis shift.
    slice_start::Int64 = 0                      # First slice index.
    slice_end::Int64 = 0                        # Last slice index.
    descrip::NTuple{80,UInt8} = ztuple(80)      # any text you like.
    aux_file::NTuple{24,UInt8} = ztuple(24)     # auxiliary filename.
    qform_code::Int32 = 0                       # NIFTI_XFORM_* code.
    sform_code::Int32 = 0                       # NIFTI_XFORM_* code.
    quatern::NTuple{3,Float32} = ztuple(3)      # Quaternion b/c/d param.
    qoffset::NTuple{3,Float32} = ztuple(3)      # Quaternion x/y/z param.
    srow_x::NTuple{4,Float64} = ztuple(4)       # 1st row affine transform.
    srow_y::NTuple{4,Float64} = ztuple(4)       # 2nd row affine transform.
    srow_z::NTuple{4,Float64} = ztuple(4)       # 3rd row affine transform.
    slice_code::Int32 = 0                       # Slice timing order.
    xyzt_units::Int32 = 0                       # Units of pixdim[1..4]
    intent_code::Int32 = 0                      # NIFTI_INTENT_* code.
    intent_name::NTuple{16,UInt8} = ztuple(16)  # 'name' or meaning of data.
    dim_info::UInt8 = 0                         # MRI slice ordering.
    unused_str::NTuple{15,UInt8} = ztuple(15)   # unused, filled with \0
end


"""
Check if a nifti_1_header struct needs to be byte swapped.
Returns 1 if it needs to be swapped, 0 if it does not.
"""
needs_swap(h::Nifti1Header) = (bswap(h.sizeof_hdr) == 348)
needs_swap(h::Nifti2Header) = (bswap(h.sizeof_hdr) == 540)


# Open stream if needed, using the correct opener
function _open_stream(path::AbstractString, mode="r")
    io = open(path, mode)
    if endswith(lowercase(path), ".gz")
        stream = mode == "r" ? GzipDecompressorStream(io) : GzipCompressorStream(io)
        return OpenedStream(stream, [stream, io])
    else
        return _open_stream(OpenedStream(io, [io]))
    end
end

function _open_stream(io::IO, mode="r")
    if mode == "r" && startswith(peek(io,String), GZIP_MAGIC)
        stream = GzipDecompressorStream(io)
        return OpenedStream(stream, [stream])
    else
        return OpenedStream(io, IO[])
    end
end

GZIP_MAGIC = 0x8b1f
GZIP_MAGIC = String(UInt8[GZIP_MAGIC >>> 8, GZIP_MAGIC & 0xff])

function _open_stream(io::OpenedStream, mode="r")
    if mode == "r" && startswith(io.stream, GZIP_MAGIC)
        stream = GzipDecompressorStream(io.stream)
        return OpenedStream(stream, [io.stream; io.mine])
    else
        return io
    end
end

_open_stream(io::OpenedStream{GzipDecompressorStream}, mode="r") = io
_open_stream(io::OpenedStream{GzipCompressorStream}, mode="w") = io


# Read header, bswap and tell if bswap needed
function _read!(io::IOType, header::NiftiHeader)
    io = _open_stream(io)
    Base.read!(io.stream, Ref(header))
    map(close, io.mine)
    do_swap = needs_swap(header)
    if do_swap
        smap!(bswap, header)
    end
    return (header=header, bswap=do_swap)
end
_read(io::IOType, type::Type{<:NiftiHeader}) = _read!(io, type())

# Read header and bswap (only return header)
read!(io::IO, header::NiftiHeader) = _read!(io, header)[:header]
read(io::IO, type::Type{<:NiftiHeader}) = _read!(io, type())[:header]
read!(io::AbstractString, header::NiftiHeader) = _read!(io, header)[:header]
read(io::AbstractString, type::Type{<:NiftiHeader}) = _read!(io, type())[:header]

# Write header
function write(io::IOType, header::NiftiHeader, bswap::Bool=false)
    io = _stream(io, 'w')
    header = bswap ? smap(Base.bswap, header) : header
    write(io.stream, Ref(header))
    map(close, io.mine)
    return io.stream
end

# Convert from NIfTI.jl to ours
function convert(::Type{Nifti1Header}, x::NIfTI.NIfTI1Header)
    ptr = reinterpret(Ptr{UInt8}, pointer_from_objref(x))
    vec = unsafe_wrap(Vector{UInt8}, ptr, sizeof(x), own=false)
    read(IOBuffer(vec), Nifti1Header)
end

# Convert from ours to NIfTI.jl
function convert(::Type{NIfTI.NIfTI1Header}, x::Nifti1Header)
    ptr = reinterpret(Ptr{UInt8}, pointer_from_objref(x))
    vec = unsafe_wrap(Vector{UInt8}, ptr, sizeof(x), own=false)
    read(IOBuffer(vec), NIfTI.NIfTI1Header)[1]
end

function bytesencode(x::NiftiHeader)
    ptr = reinterpret(Ptr{UInt8}, pointer_from_objref(x))
    vec = unsafe_wrap(Vector{UInt8}, ptr, sizeof(x), own=false)
end

function bytesdecode(::Type{Nifti1Header}, x::Vector{UInt8})
    read(IOBuffer(x), Nifti1Header)
end

function bytesdecode(::Type{NIfTI.NIfTI1Header}, x::Vector{UInt8})
    read(IOBuffer(x), NIfTI.NIfTI1Header)[1]
end

function base64encode(x::NiftiHeader)
    ptr = reinterpret(Ptr{UInt8}, pointer_from_objref(x))
    vec = unsafe_wrap(Vector{UInt8}, ptr, sizeof(x), own=false)
    base64encode(vec)
end

function base64decode(x::AbstractString, T::Type{<:NiftiHeader})
    read(Base64DecodePipe(IOBuffer(x)), T)
end

@enumx DataTypeCode begin
    u1 = 2
    i2 = 4
    i4 = 8
    f4 = 16
    c8 = 32
    rgb = 128
    i1 = 256
    u2 = 512
    u4 = 768
    i8 = 1024
    u8 = 1280
    f16 = 1536
    c16 = 1792
    c32 = 2048
    rgba = 2304
end

@enumx UnitCode begin
    unknown = 0
    meter = 1
    millimeter = 2
    micrometer = 3
    second = 8
    millisecond = 16
    microsecond = 24
    hertz = 32
    micro = 40
    radian = 48
end

@enumx IntentCode begin
    NONE = 0
    CORREL = 2
    TTEST = 3
    FTEST = 4
    ZSCORE = 5
    CHISQ = 6
    BETA = 7
    BINOM = 8
    GAMMA = 9
    POISSON = 10
    NORMAL = 11
    FTEST_NONC = 12
    CHISQ_NONC = 13
    LOGISTIC = 14
    LAPLACE = 15
    UNIFORM = 16
    TTEST_NONC = 17
    WEIBULL = 18
    CHI = 19
    INVGAUSS = 20
    EXTVAL = 21
    PVAL = 22
    LOGPVAL = 23
    LOG10PVAL = 24
    ESTIMATE = 1001
    LABEL = 1002
    NEURONAME = 1003
    GENMATRIX = 1004
    SYMMATRIX = 1005
    DISPVECT = 1006
    VECTOR = 1007
    POINTSET = 1008
    TRIANGLE = 1009
    QUATERNION = 1010
    DIMLESS = 1011
    TIME_SERIES = 2001
    NODE_INDEX = 2002
    RGB_VECTOR = 2003
    RGBA_VECTOR = 2004
    SHAPE = 2005
    FSL_FNIRT_DISPLACEMENT_FIELD = 2006
    FSL_CUBIC_SPLINE_COEFFICIENTS = 2007
    FSL_DCT_COEFFICIENTS = 2008
    FSL_QUADRATIC_SPLINE_COEFFICIENTS = 2009
    FSL_TOPUP_CUBIC_SPLINE_COEFFICIENTS = 2016
    FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS =  2017
    FSL_TOPUP_FIELD = 2018
end

IntentNbPrm = Dict(key => 0 for key in instances(IntentCode.T))
merge!(IntentNbPrm, Dict(
    IntentCode.CORREL => 1,
    IntentCode.TTEST => 1,
    IntentCode.FTEST => 2,
    IntentCode.CHISQ => 1,
    IntentCode.BETA => 2,
    IntentCode.BINOM => 2,
    IntentCode.GAMMA => 2,
    IntentCode.POISSON => 1,
    IntentCode.NORMAL => 2,
    IntentCode.FTEST_NONC => 3,
    IntentCode.CHISQ_NONC => 2,
    IntentCode.LOGISTIC => 2,
    IntentCode.LAPLACE => 2,
    IntentCode.UNIFORM => 2,
    IntentCode.TTEST_NONC => 2,
    IntentCode.WEIBULL => 2,
    IntentCode.CHI => 1,
    IntentCode.INVGAUSS => 2,
    IntentCode.EXTVAL => 2,
    IntentCode.GENMATRIX => 2,
    IntentCode.SYMMATRIX => 1,
))

@enumx XFormCode begin
    UNKNOWN = 0
    SCANNER_ANAT = 1
    ALIGNED_ANAT = 2
    TALAIRACH = 3
    MNI = 4
    TEMPLATE_OTHER = 5
end

@enumx SliceOrderCode begin
    UNKNOWN = 0
    SEQ_INC = 1
    SEQ_DEC = 2
    ALT_INC = 3
    ALT_DEC = 4
    ALT_INC2 = 5
    ALT_DEC2 = 6
end

UnitRecoder(x::Integer) = string(UnitCode.T(x))
UnitRecoder(x::String) = getproperty(UnitCode, Symbol(x))
IntentRecoder(x::Integer) = string(IntentCode.T(x))
IntentRecoder(x::String) = getproperty(IntentCode, Symbol(x))
XFormRecoder(x::Integer) = string(XFormCode.T(x))
XFormRecoder(x::String) = getproperty(XFormCode, Symbol(x))
SliceOrderRecoder(x::Integer) = string(SliceOrderCode.T(x))
SliceOrderRecoder(x::String) = getproperty(SliceOrderCode, Symbol(x))

function DataTypeRecoder(x::Integer)
    if x == DataTypeCode.rgb
        (("r", "u1"), ("g", "u1"), ("b", "u1"))
    elseif x == DataTypeCode.rgba
        (("r", "u1"), ("g", "u1"), ("b", "u1"), ("a", "u1"))
    else
        string(DataTypeCode.T(x))
    end
end
DataTypeRecoder(x::String) = getproperty(DataTypeCode, Symbol(x))
function DataTypeRecoder(x::Tuple)
    if x == (("r", "u1"), ("g", "u1"), ("b", "u1"))
        DataTypeRecoder("rgb")
    elseif x == (("r", "u1"), ("g", "u1"), ("b", "u1"), ("a", "u1"))
        DataTypeRecoder("rgba")
    else
        error("unsupported data type")
    end
end
function DataTypeRecoder(x::Array)
    DataTypeRecoder(map(Tuple, Tuple(x)))
end
