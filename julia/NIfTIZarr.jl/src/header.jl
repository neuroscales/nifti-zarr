import Base: @kwdef, read, write
import EnumX: @enumx
import Printf: @printf
import Base64: base64encode, Base64DecodePipe
import NIfTI: NIfTIHeader, NIfTI1Header, NIfTI2Header
# import GZip
import CodecZlib: GzipCompressorStream, GzipDecompressorStream
import Bijections: Bijection

"""A generic IO type"""
IOType = Union{IO,AbstractString}

struct OpenedStream{T<:IO}
    stream::T
    mine::Vector{<:IO}
end



function Base.show(x::NIfTIHeader)
    fieldlength = maximum(map(f -> length(string(f)), fieldnames(typeof(x)))) + 1
    print("$(typeof(x))()\n")
    for field in fieldnames(typeof(x))
        @printf "%s => %s\n" field reprfield(getproperty(x, field))
    end
end

"""
    ztuple(n::Integer)

Create a NTuple of zeros
"""
ztuple(n) = ntuple(x -> 0, n)

# Valid magic strings
MAGIC1 = UInt8.(('n', 'i', '1', '\0'))   # two files
MAGIC1P = UInt8.(('n', '+', '1', '\0'))  # single file
MAGIC2 = UInt8.(('n', 'i', '2', '\0', '\u0D', '\u0A', '\u1D', '\u1A'))   # two files
MAGIC2P = UInt8.(('n', '+', '2', '\0', '\u0D', '\u0A', '\u1D', '\u1A'))  # single file

"""
Check if a nifti_1_header struct needs to be byte swapped.
Returns 1 if it needs to be swapped, 0 if it does not.
"""
needs_swap(h::NIfTI1Header) = (bswap(h.sizeof_hdr) == 348)
needs_swap(h::NIfTI2Header) = (bswap(h.sizeof_hdr) == 540)


# Open stream if needed, using the correct opener
function _open_stream(path::AbstractString, mode="r")
    io = open(path, mode)
    if endswith(lowercase(path), ".gz")
        stream = mode == "r" ? GzipDecompressorStream(io) : GzipCompressorStream(io)
        return OpenedStream(stream, [stream, io])
    else
        return OpenedStream(io, [io])
    end
end

function _open_stream(io::IO, mode="r")
    if mode == "r" && startswith(peek(io, String), GZIP_MAGIC)
        stream = GzipDecompressorStream(io)
        return OpenedStream(stream, [stream])
    else
        return OpenedStream(io, IO[])
    end
end

GZIP_MAGIC = 0x8b1f
GZIP_MAGIC = String(UInt8[GZIP_MAGIC>>>8, GZIP_MAGIC&0xff])

function _open_stream(io::OpenedStream, mode="r")
    if mode == "r" && startswith(peek(io.stream), GZIP_MAGIC)
        stream = GzipDecompressorStream(io.stream)
        return OpenedStream(stream, [io.stream; io.mine])
    else
        return io
    end
end

_open_stream(io::OpenedStream{GzipDecompressorStream}, mode="r") = io
_open_stream(io::OpenedStream{GzipCompressorStream}, mode="w") = io


# Read header, bswap and tell if bswap needed
function _read!(io::IOType, header::NIfTIHeader)
    io = _open_stream(io)
    Base.read!(io.stream, Ref(header))
    map(close, io.mine)
    do_swap = needs_swap(header)
    if do_swap
        smap!(bswap, header)
    end
    return (header=header, bswap=do_swap)
end
_read(io::IOType, type::Type{<:NIfTIHeader}) = _read!(io, type())

# Read header and bswap (only return header)
read!(io::IO, header::NIfTIHeader) = _read!(io, header)[:header]
read(io::IO, type::Type{<:NIfTIHeader}) = _read!(io, type())[:header]
read!(io::AbstractString, header::NIfTIHeader) = _read!(io, header)[:header]
read(io::AbstractString, type::Type{<:NIfTIHeader}) = _read!(io, type())[:header]

# Write header
function write(io::IOType, header::NIfTIHeader, bswap::Bool=false)
    io = _stream(io, 'w')
    header = bswap ? smap(Base.bswap, header) : header
    write(io.stream, Ref(header))
    map(close, io.mine)
    return io.stream
end

function bytesencode(x::NIfTIHeader)
    ptr = reinterpret(Ptr{UInt8}, pointer_from_objref(x))
    vec = unsafe_wrap(Vector{UInt8}, ptr, sizeof(x), own=false)
end

function bytesdecode(::Type{NIfTI1Header}, x::Vector{UInt8})
    read(IOBuffer(x), Nifti1Header)
end

function bytesdecode(::Type{NIfTI1Header}, x::Vector{UInt8})
    read(IOBuffer(x), NIfTI1Header)[1]
end

function base64encode(x::NIfTIHeader)
    ptr = reinterpret(Ptr{UInt8}, pointer_from_objref(x))
    vec = unsafe_wrap(Vector{UInt8}, ptr, sizeof(x), own=false)
    base64encode(vec)
end

function base64decode(x::AbstractString, T::Type{<:NIfTIHeader})
    read(Base64DecodePipe(IOBuffer(x)), T)
end

DataTypeCode = Dict(
    "uint8" => 2,
    "int16" => 4,
    "int32" => 8,
    "single" => 16,
    "complex64" => 32,
    "double" => 64,
    "rgb" => 128,
    "int8" => 256,
    "uint16" => 512,
    "uint32" => 768,
    "int64" => 1024,
    "uint64" => 1280,
    "double128" => 1536,
    "complex128" => 1792,
    "complex256" => 2048,
    "rgba" => 2304,
)

UnitCode = Dict(
    "" => 0,
    "m" => 1,
    "mm" => 2,
    "um" => 3,
    "s" => 8,
    "ms" => 16,
    "us" => 24,
    "hz" => 32,
    "ppm" => 40,
    "rad/s" => 48
)

JNIfTIUnit2ZarrUnit = Dict(
    "" => "",
    "m" => "meter",
    "mm" => "millimeter",
    "um" => "micrometer",
    "s" => "second",
    "ms" => "millisecond",
    "us" => "microsecond",
    "hz" => "hertz",
    "ppm" => "micro",
    "rad/s" => "radian"
)

JNIfTIDataType2ZarrDataType = Dict(
    "uint8" => "u1",
    "int16" => "i2",
    "int32" => "i4",
    "single" => "f4",
    "complex64" => "c8",
    (("r", "uint8"), ("g", "uint8"), ("b", "uint8")) => (("r", "u1"), ("g", "u1"), ("b", "u1")),
    "int8" => "i1",
    "uint16" => "u2",
    "uint32" => "u4",
    "int64" => "i8",
    "uint64" => "u8",
    "double128" => "f16",
    "complex128" => "c16",
    "complex256" => "c32",
    (("r", "uint8"), ("g", "uint8"), ("b", "uint8"), ("a", "uint8")) => (("r", "u1"), ("g", "u1"), ("b", "u1"), ("a", "u1")),
)



IntentCode = Dict(
    "" => 0,
    "corr" => 2,
    "ttest" => 3,
    "ftest" => 4,
    "zscore" => 5,
    "chi2" => 6,
    "beta" => 7,
    "binomial" => 8,
    "gamma" => 9,
    "poisson" => 10,
    "normal" => 11,
    "ncftest" => 12,
    "ncchi2" => 13,
    "logistic" => 14,
    "laplace" => 15,
    "uniform" => 16,
    "ncttest" => 17,
    "weibull" => 18,
    "chi" => 19,
    "invgauss" => 20,
    "extval" => 21,
    "pvalue" => 22,
    "logpvalue" => 23,
    "log10pvalue" => 24,
    "estimate" => 1001,
    "label" => 1002,
    "neuronames" => 1003,
    "matrix" => 1004,
    "symmatrix" => 1005,
    "dispvec" => 1006,
    "vector" => 1007,
    "point" => 1008,
    "triangle" => 1009,
    "quaternion" => 1010,
    "unitless" => 1011,
    "tseries" => 2001,
    "elem" => 2002,
    "rgb" => 2003,
    "rgba" => 2004,
    "shape" => 2005,
    "fsl_fnirt_displacement_field" => 2006,
    "fsl_cubic_spline_coefficients" => 2007,
    "fsl_dct_coefficients" => 2008,
    "fsl_quadratic_spline_coefficients" => 2009,
    "fsl_topup_cubic_spline_coefficients" => 2016,
    "fsl_topup_quadratic_spline_coefficients" => 2017,
    "fsl_topup_field" => 2018
)


IntentNbPrm = Dict(key => 0 for (key, v) in IntentCode)
merge!(IntentNbPrm, Dict(
    "corr" => 1,
    "ttest" => 1,
    "ftest" => 2,
    "chi2" => 1,
    "beta" => 2,
    "binomial" => 2,
    "gamma" => 2,
    "poisson" => 1,
    "normal" => 2,
    "ncftest" => 3,
    "ncchi2" => 2,
    "logistic" => 2,
    "laplace" => 2,
    "uniform" => 2,
    "ncttest" => 2,
    "weibull" => 2,
    "chi" => 1,
    "invgauss" => 2,
    "extval" => 2,
    "matrix" => 2,
    "symmatrix" => 1,
))

XFormCode = Dict(
    "" => 0,
    "scanner_anat" => 1,
    "aligned_anat" => 2,
    "talairach" => 3,
    "mni_152" => 4,
    "template_other" => 5
)

SliceOrderCode = Dict(
    "" => 0,
    "seq+" => 1,
    "seq-" => 2,
    "alt+" => 3,
    "alt-" => 4,
    "alt2+" => 5,
    "alt2-" => 6
)


DataTypeCode = Bijection(DataTypeCode)
UnitCode = Bijection(UnitCode)
IntentCode = Bijection(IntentCode)
XFormCode = Bijection(XFormCode)
SliceOrderCode = Bijection(SliceOrderCode)

UnitRecoder(x::Integer) = UnitCode(x)
UnitRecoder(x::String) = UnitCode[x]
IntentRecoder(x::Integer) = IntentCode(x)
IntentRecoder(x::String) = IntentCode[x]
XFormRecoder(x::Integer) = XFormCode(x)
XFormRecoder(x::String) = XFormCode[x]
SliceOrderRecoder(x::Integer) = SliceOrderCode(x)
SliceOrderRecoder(x::String) = SliceOrderCode[x]

function DataTypeRecoder(x::Integer)
    if x == DataTypeCode["rgb"]
        (("r", "uint8"), ("g", "uint8"), ("b", "uint8"))
    elseif x == DataTypeCode["rgba"]
        (("r", "uint8"), ("g", "uint8"), ("b", "uint8"), ("a", "uint8"))
    else
        DataTypeCode(x)
    end
end
DataTypeRecoder(x::String) = DataTypeCode[x]
function DataTypeRecoder(x::Tuple)
    if x == (("r", "uint8"), ("g", "uint8"), ("b", "uint8"))
        DataTypeRecoder["rgb"]
    elseif x == (("r", "uint8"), ("g", "uint8"), ("b", "uint8"), ("a", "uint8"))
        DataTypeRecoder["rgba"]
    else
        error("unsupported data type")
    end
end
function DataTypeRecoder(x::Array)
    DataTypeRecoder(map(Tuple, Tuple(x)))
end
