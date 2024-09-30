import EnumX: @enumx
import Images: gaussian_pyramid
import NIfTI
import Zarr
import Base: delete!

@enumx CompressorType no = raw = 0 blosc zlib
@enumx PyramidMethod gaussian
@enumx IsLabel yes no auto
@enumx ByteOrder little big
SystemByteOrder = ENDIAN_BOM == 0x04030201 ? ByteOrder.little : ByteOrder.big
ByteOrderSymbol = Dict(ByteOrder.little => '<', ByteOrder.big => '>')

_tuple2string(x::NTuple) = split(join(Char.(x)), '\0')[begin]

nii2json(header::NiftiHeader) = nii2json(header, SystemByteOrder)

"""
    nii2json(header[, byteorder])

Convert a nifti header into a JSON-compatible dictionary

# Example
```julia-repl
julia> import NIfTIZarr: NiftiHeader, read, nii2json
julia> header = read("path/to/nifti.nii.gz", NiftiHeader)
julia> jsonheader = nii2json(header)
```
"""
function nii2json(header::NiftiHeader, byteorder::ByteOrder.T)

    magic = _tuple2string(header.magic)
    magic = magic[1] * 'z' * magic[3]
    ndim = header.dim[1]
    intent_code = IntentRecoder(header.intent_code)
    nb_intent_prm = IntentNbPrm[IntentCode.T(header.intent_code)]

    jsonheader = Dict(
        "magic" => magic,
        "dim" => collect(header.dim[2:2+ndim-1]),
        "pixdim" => collect(header.pixdim[2:2+ndim-1]),
        "units" => Dict(
            "space" => UnitRecoder(header.xyzt_units & 0x07),
            "time" => UnitRecoder(header.xyzt_units & 0x38)
        ),
        "datatype" => DataTypeRecoder(header.datatype),
        "dim_info" => Dict(
            "freq" => header.dim_info & 0x03,
            "phase" => (header.dim_info >> 2) & 0x03,
            "slice" => (header.dim_info >> 4) & 0x03
        ),
        "intent" => Dict(
            "code" => intent_code,
            "name" => _tuple2string(header.intent_name),
            "p" => collect(header.intent_p[begin:nb_intent_prm])
        ),
        "scl" => Dict(
            "slope" => header.scl_slope,
            "inter" => header.scl_inter
        ),
        "slice" => Dict(
            "code" => SliceOrderRecoder(header.slice_code),
            "start" => header.slice_start,
            "end" => header.slice_end,
            "duration" => header.slice_duration
        ),
        "cal" => Dict(
            "min" => header.cal_min,
            "max" => header.cal_max
        ),
        "toffset" => header.toffset,
        "descrip" => _tuple2string(header.descrip),
        "aux_file" => _tuple2string(header.aux_file),
        "qform" => Dict(
            "intent" => XFormRecoder(header.qform_code),
            "quatern" => collect(header.quatern),
            "offset" => collect(header.qoffset),
            "fac" => header.pixdim[begin]
        ),
        "sform" => Dict(
            "intent" => XFormRecoder(header.sform_code),
            "affine" => [[header.srow_x...] [header.srow_y...] [header.srow_z...]]
        )
    )
    if !isfinite(jsonheader["scl"]["slope"])
        jsonheader["scl"]["slope"] = 0.0
    end
    if !isfinite(jsonheader["scl"]["inter"])
        jsonheader["scl"]["inter"] = 0.0
    end

    # Fix data type
    bo = ByteOrderSymbol[byteorder]
    if isa(jsonheader["datatype"], Array)
        jsonheader["datatype"] = map(
            x -> [x[1], '|' * x[2]],
            jsonheader["datatype"]
        )
    elseif jsonheader["datatype"][end] == '1'
        jsonheader["datatype"] = '|' + jsonheader["datatype"]
    else
        jsonheader["datatype"] = bo * jsonheader["datatype"]
    end

    return jsonheader
end


function _make_pyramid3d(data3d::AbstractArray, nb_levels::Integer, label::Bool=false)
    nbatch = size(data3d)[1:end-3]
    nxyz = size(data3d)[end-2:end] 
    data4d = reshape(data3d, (prod(nbatch), nxyz...))
    if nb_levels < 0
        T = typeof(nb_levels)
        nb_levels = min(T.(floor.(log2.(nxyz) .- 2))...)
    end
    max_layer = nb_levels - 1

    # I would really prefer to use sigma = 0.42, which corresponds to a
    # full-width at half-maximum of 1 (which, assuming that the signal at
    # the previous resolution has a PSF of 1 pixel, ensures that the filtered
    # image has a PSF of 2 voxels, matching the downsampling rate).
    # However, gaussian_pyramid uses a IIF filter and complains when sigma < 1
    pyramid_values(x::AbstractArray) = gaussian_pyramid(x, max_layer, 2, 1.0)

    function pyramid_labels(x::AbstractArray)
        labels = unique(x)
        pyramid_maxprob = pyramid_values(x == labels[begin])
        pyramid = map(y -> fill!(similar(y, eltype(x)), 0), pyramid_maxprob)

        for label in labels[2:end]
            pyramid_prob = pyramid_values(x == label)
            for (value, prob, maxprob) in zip(pyramid, pyramid_prob, pyramid_maxprob)
                mask = prob > maxprob
                value[mask] = label[mask]
                maxprob[mask] = prob[mask]
            end
        end
        return pyramid
    end

    pyramids3d = map(label ? pyramid_labels : pyramid_values, eachslice(data4d; dims=1))

    level0 = map(x->reshape(x[begin], (1, size(x[begin])...)), pyramids3d)
    level0 = reduce(vcat, level0)
    level0 = reshape(level0, (nbatch..., size(level0)[2:end]...))
    levels = [level0]

    for n in 2:length(pyramids3d[1])
        level = reduce(vcat, map(x->reshape(x[n], (1, size(x[n])...)), pyramids3d))
        level = reshape(level, (nbatch..., size(level)[2:end]...))
        push!(levels, level)
    end

    return levels
end

_compressor_map = Dict(
    CompressorType.no => Zarr.NoCompressor,
    CompressorType.blosc => Zarr.BloscCompressor,
    CompressorType.zlib => Zarr.ZlibCompressor,
)
_make_compressor(x; kwargs...) = x
_make_compressor(name::CompressorType.T; kwargs...) = _compressor_map[name](;kwargs...)

_nii2zarr_prepare_inp(x) = NIfTI.niread(x)
_nii2zarr_prepare_inp(x::NIfTI.NIVolume) = x
function _nii2zarr_prepare_out(x)
    try
        Zarr.zgroup(x)
    catch
        Zarr.zopen(x, "w")
    end
end
_nii2zarr_prepare_out(x::Zarr.ZGroup) = x
_nii2zarr_prepare_chunk(x::NTuple{N, <:Integer} where N) = (x[begin:3],)
_nii2zarr_prepare_chunk(x::Integer) = [(x, x, x, 1, 1)]
_nii2zarr_prepare_chunk(x::NTuple{1, <:Integer}) = [(x..., x..., x..., 1, 1)]
_nii2zarr_prepare_chunk(x::NTuple{2, <:Integer}) = [(x..., x[end], 1, 1)]
_nii2zarr_prepare_chunk(x::NTuple{3, <:Integer}) = [(x, 1, 1)]
_nii2zarr_prepare_chunk(x::NTuple{4, <:Integer}) = [(x, 1)]
_nii2zarr_prepare_chunk(x::NTuple{5, <:Integer}) = [x]
_nii2zarr_prepare_chunk(x::Tuple{Tuple}) = _nii2zarr_prepare_chunk(collect(x))
_nii2zarr_prepare_chunk(x::Vector) = map(_nii2zarr_prepare_chunk, x)
_nii2zarr_prepare_chunk(x::Vector{NTuple{5, T}} where T <: Integer) = x
_nii2zarr_prepare_label(x::String) = getfield(IsLabel, Symbol(x))
_nii2zarr_prepare_label(x::Symbol) = getfield(IsLabel, x)
_nii2zarr_prepare_label(x::IsLabel.T) = x
_nii2zarr_prepare_compressor(x::String) = getfield(CompressorType, Symbol(x))
_nii2zarr_prepare_compressor(x::Symbol) = getfield(CompressorType, x)
_nii2zarr_prepare_compressor(x::CompressorType.T) = x
_nii2zarr_prepare_method(x::String) = getfield(MethodType, Symbol(x))
_nii2zarr_prepare_method(x::Symbol) = getfield(MethodType, x)
_nii2zarr_prepare_method(x::PyramidMethod.T) = x


function _nii2zarr_prepare_chunk(x::Array, nb_levels)
    if length(x) > nb_levels
        return _nii2zarr_prepare_chunk(x[1:nb_levels])
    else
        while length(x) < nb_levels
            push!(x, x[end])
        end
        return _nii2zarr_prepare_chunk(x)
    end
end


function _nii2zarr_ome_attrs(axes, shapes, niiheader)
    attrs = Dict(
        "multiscales" => [Dict()]
    )
    multiscales = attrs["multiscales"][1]

    typemap = Dict(
        "t" => "time",
        "c" => "channel",
        "z" => "space",
        "y" => "space",
        "x" => "space",
    )
    idxmap = Dict(
        "x" => 1,
        "y" => 2,
        "z" => 3,
        "t" => 4,
    )

    multiscales["axes"] = [(
        typemap[axis] == "channel"
        ? Dict(
            "name" => axis,
            "type" => typemap[axis]
        )
        : Dict(
            "name" => axis,
            "type" => typemap[axis],
            "unit" => niiheader["units"][typemap[axis]]
        )
    ) for axis in axes]

    multiscales["datasets"] = [Dict() for _ in 1:length(shapes)]
    for n in eachindex(shapes)
        # spatial scale at each pyramid level
        level = multiscales["datasets"][n]

        # Image.gaussian_pyramid aligns voxel edges across scales
        # (same behavior as scipy.ndimage.zoom(..., grid_mode=True))
        # so the effective scaling is the shape ratio, and there is
        # a half voxel shift wrt to the "center of first voxel" frame
        level["coordinateTransformations"] = [
            Dict(
                "type" => "scale",
                "scale" => [(
                    typemap[axis] == "spatial"
                    ? (shapes[1][i]/shapes[n][i])*niiheader["pixdim"][idxmap[axis]]
                    : 1.0
                ) for (i, axis) in enumerate(axes)]
            ),
            Dict(
                "type" => "translation",
                "translation" => [(
                    typemap[axis] == "spatial"
                    ? (shapes[1][i]/shapes[n][i] - 1)*niiheader["pixdim"][idxmap[axis]]*0.5
                    : 0.0
                ) for (i, axis) in enumerate(axes)]
            ),
        ]
    end

    # time scale
    multiscales["coordinateTransformations"] = [Dict(
        "scale" => [(
            typemap[axis] == "time"
            ? niiheader["pixdim"][idxmap[axis]]
            : 1.0
        )  for (i, axis) in enumerate(axes)],
        "type" => "scale"
    )]

    return attrs
end


# massage inputs that have a "liberal" format
nii2zarr(inp, out;
    chunk = 64,
    method = PyramidMethod.gaussian,
    label = IsLabel.auto,
    compressor = CompressorType.blosc,
    other_options...
) = nii2zarr(
    _nii2zarr_prepare_inp(inp),
    _nii2zarr_prepare_out(out);
    chunk = _nii2zarr_prepare_chunk(chunk),
    method = _nii2zarr_prepare_method(method),
    label = _nii2zarr_prepare_label(label),
    compressor = _nii2zarr_prepare_compressor(compressor),
    other_options...
)


"""
    nii2zarr(inp, out; [chunk, nb_levels, method, label, fill_value, compressor])

Convert a nifti file to nifti-zarr directory.

- The input should be a `NIfTI.NIVolume`, but can also be anything that
  `NIfTI.niread` handles.
- The output should be a `Zarr.ZGroup`, but can also be anything that
  `Zarr.zgroup` handles.

# Keyword arguments
- `chunk`: Chunk size, per x/y/z/t/c, per level.
    * The inner tuple allows different chunk sizes to be used along
      each dimension.
    * The outer array allows different chunk sizes to be used at
      different pyramid levels.
    * By default, the chunk size for axes "t" and "c" is one.
- `nb_levels`: Number of levels in the pyramid. Default: all possible levels.
- `method`: Method used to compute the pyramid. Only gaussian supported.
- `label`: Is this is a label volume?  By default, guess from intent code.
- `fill_value`: Value to use for missing tiles
- `compressor`: Compression to use. Only "blosc" and "raw" implemented.

# Example
```julia-repl
julia> import NIfTIZarr: nii2zarr
julia> nii2zarr("path/to/nifti.nii.gz", "s3://path/to/bucket")
```
"""
function nii2zarr(inp::NIfTI.NIVolume, out::Zarr.ZGroup;
    chunk::Vector{NTuple{5, T}} where T <: Integer = [(64, 64, 64, 1, 1)],
    nb_levels::Integer = -1,
    method::PyramidMethod.T = PyramidMethod.gaussian,
    label::IsLabel.T = IsLabel.auto,
    fill_value = nothing,
    compressor::CompressorType.T = CompressorType.blosc,
    compressor_options...
)
    # Convert NIfTI.jl's header into my header (they only handle nifti-1)
    header = convert(Nifti1Header, inp.header)
    jsonheader = nii2json(header)

    # Fix array shape
    if ndims(inp.raw) == 3
        nbatch = 0
        perm = [3, 2, 1]
        axes = ["z", "y", "x"]
        ARRAY_DIMENSIONS = ["z", "y", "x"]
    elseif ndims(inp.raw) == 4
        nbatch = 1
        perm = [4, 3, 2, 1]
        axes = ["t", "z", "y", "x"]
        ARRAY_DIMENSIONS = ["time", "z", "y", "x"]
    elseif ndims(inp.raw) == 5
        nbatch = 2
        perm = [4, 5, 3, 2, 1]
        axes = ["t", "c", "z", "y", "x"]
        ARRAY_DIMENSIONS = ["time", "channel", "z", "y", "x"]
    elseif ndims(inp.raw) > 5
        error("Too many dimensions for conversion to nii.zarr")
    else
        error("Too few dimensions for conversion to nii.zarr. Is this really an image?")
    end
    data = PermutedDimsArray(inp.raw, perm)

    # Compute image pyramid
    bool_label = jsonheader["intent"]["code"] in ("LABEL", "NEURONAMES")
    if label == IsLabel.yes
        bool_label = true
    elseif label == IsLabel.no
        bool_label = false
    end

    data = _make_pyramid3d(data, nb_levels, bool_label)
    nb_levels = length(data)
    shapes = [size(d) for d in data]
    attrs = _nii2zarr_ome_attrs(axes, shapes, jsonheader)

    # Prepare array metadata at each level
    compressor = _make_compressor(compressor; compressor_options...)
    chunk = _nii2zarr_prepare_chunk(chunk, nb_levels)
    chunk = map(x->Dict("x"=>x[1], "y"=>x[2], "z"=>x[3], "t"=>x[4], "c"=>x[5]), chunk)
    chunk = map(x->Tuple(x[axis] for axis in axes), chunk)

    # Write nifti header
    header.magic = MAGIC1Z
    binheader = bytesencode(header)
    delete!(out.storage, "", "nifti")
    header_array = Zarr.zcreate(
        eltype(binheader), out, "nifti", length(binheader);
        chunks=(length(binheader),),
        fill_as_missing=false,
        fill_value=nothing,
        compressor=Zarr.NoCompressor(),
        attrs=jsonheader
    )
    # header_array[:] = binheader
    copy!(header_array, binheader)

    # Write group attributes
    Zarr.writeattrs(out.storage, out.path, attrs)

    # Write zarr arrays
    for n in eachindex(shapes)
        delete!(out.storage, "", string(n-1))
        subarray = Zarr.zcreate(
            eltype(data[n]), out, string(n-1), shapes[n]...;
            chunks=chunk[n],
            fill_as_missing=false,
            fill_value=nothing,
            compressor=compressor,
            # dimension_separator='/',
            # order='F',
            attrs = Dict("_ARRAY_DIMENSIONS" => ARRAY_DIMENSIONS)
        )
        # subarray[:,:,:] = data[n]
        copy!(subarray, data[n])
    end

end


# Fix delete! from Zarr.jl, which implements it as:
#   Base.delete!(s::DirectoryStore, k::String) = isfile(joinpath(s.folder, k)) && rm(joinpath(s.folder, k))
# Their version only deletes files, and not directories.
delete!(s::Zarr.DirectoryStore, k::String) = rm(joinpath(s.folder, k); force=true, recursive=true)
