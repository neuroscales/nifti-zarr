import EnumX: @xenum
import Base64: base64encode
import Images: gaussian_pyramid
import NIfTI
import Zarr

@xenum CompressorType no = raw = 0 blosc zlib
@xenum PyramidMethod gaussian
@xenum IsLabel yes no auto
@xenum ByteOrder little big
SystemByteOrder = (ENDIAN_BOM == 0x04030201) ? ByteOrder.little : ByteOrder.big
ByteOrderSymbol = Dict(ByteOrder.little => '<', ByteOrder.big => '>')


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
function nii2json(header::NiftiHeader, byteorder::ByteOrder = SystemByteOrder)

    magic = String(header.magic)
    magic = magic[1:1] * "z" * magic[2:end]
    ndim = header.dim[1]
    intent_code = IntentRecoder(header.intentcode)
    nb_intent_prm = IntentNbPrm[header.intentcode]

    jsonheader = Dict(
        "magic" => magic,
        "dim" => collect(header.dim[2:2+ndim-1]),
        "pixdim" => collect(header.pixdim[2:2+ndim]-1),
        "units" => Dict(
            "space" => UnitRecoder(header.xyzt_units & 0x07),
            "time" => UnitRecoder(header.xyzt_units & 0x38),
        ),
        "datatype" => DataTypeRecoder(header.datatype),
        "dim_info" => Dict(
            "freq" => header.dim_info & 0x03,
            "phase" => (header.dim_info >> 2) & 0x03,
            "slice" => (header.dim_info >> 4) & 0x03,
        ),
        "intent" => Dict(
            "code" => intent_code,
            "name" => String(header.intent_name),
            "p" => collect(header.intent_p[1:1+nb_intent_prm]),
        ),
        "scl" => Dict(
            "slope" => header.scl_slope,
            "inter" => header.scl_inter,
        ),
        "slice": Dict(
            "code" => SliceOrderRecoder(header.slice_code),
            "start" => header.slice_start,
            "end" => header.slice_end,
            "duration" => header.slice_duration,
        ),
        "cal" => Dict(
            "min" => header.cal_min,
            "max" => header.cal_max,
        ),
        "toffset" => header.toffset,
        "descrip" => String(header.descrip),
        "aux_file" => String(header.aux_file),
        "qform" => Dict(
            "intent": XFormRecoder(header.qform_code),
            "quatern": collect(header.quatern),
            "offset": collect(header.qoffset),
            "fac": header.pixdim[begin],
        ),
        "sform" => Dict(
            "intent": XFormRecoder(header.sform_code),
            "affine": [header.srow_x... ; header.srow_y... ; header.srow_z...],
        ),
    )
    if !isfinite(jsonheader["scl"]["slope"])
        jsonheader["scl"]["slope"] = 0.0
    end
    if !isfinite(jsonheader["scl"]["inter"])
        jsonheader["scl"]["inter"] = 0.0
    end

    # Fix data type
    bo = ByteOrderSymbol[byteorder]
    if jsonheader["datatype"] <: Array
        jsonheader["datatype"] = map(
            x -> [x[1], '|' * x[2]],
            jsonheader["datatype"]
        )
    elseif jsonheader["datatype"][end] == '1'
        jsonheader["datatype"] = '|' + jsonheader["datatype"]
    else
        jsonheader["datatype"] = bo * jsonheader["datatype"]
    end

    # Dump the binary header
    jsonheader["base64"] = String(base64encode(header))

    return jsonheader
end


function _make_pyramid5d(data5d::AbstractArray, nb_levels::Integer, label::Bool=false)
    (nt, nc, nxyz...) = size(data5d)
    data4d = resize(data5d, (nt*nc, nxyz...))
    max_layer = nb_levels > 0 ? nb_levels - 1 : -1

    pyramid_values(x::AbstractArray) = gaussian_pyramid(x, max_layer, 2)

    function pyramid_labels(x::AbstractArray)
        labels = unique(x)
        pyramid_maxprob = gaussian_pyramid(x == labels[0], max_layer, 2)
        pyramid = map(y -> fill!(similar(y, eltype(x)), 0), pyramid_maxprob)

        for label in labels[2:end]
            pyramid_prob = gaussian_pyramid(x == label, max_layer, 2)
            for (value, prob, maxprob) in zip(pyramid, pyramid_prob, pyramid_maxprob)
                mask = prob > maxprob
                value[mask] = label[mask]
                maxprob[mask] = prob[mask]
            end
        end
        return pyramid
    end

    pyramids3d = map(label ? pyramid_labels : pyramid_values, data4d)

    level0 = stack(map(x->x[begin], pyramids3d))
    level0 = resize(level0, (nt, nc, size(level0)...))
    levels = [level0]

    for n in 2:len(pyramids3s[1])
        level = stack(map(x->x[n], pyramids3d))
        level = resize(level, (nt, nc, size(level)...))
        push!(levels, level)
    end

    return levels
end


function _make_compressor(x; kwargs...)
    return x
end


function _make_compressor(name::CompressorType; kwargs...)
    name = lower(name)
    if name == CompressorType.no
        Compressor = Zarr.NoCompressor
    elseif name == CompressorType.blosc
        Compressor = numcodecs.BloscCompressor
    elseif name == "zlib"
        error("Zlib compression not implemented")
    else
        error("Unknown compressor $name")
    end
    return Compressor(;kwargs...)
end


_nii2zarr_prepare_inp(x) = Nifti.read(x)
_nii2zarr_prepare_inp(x::NIfTI.NIVolume) = x
_nii2zarr_prepare_out(x) = Zarr.zopen(x)
_nii2zarr_prepare_out(x::Zarr.AbstractStore) = Zarr.zgroup(x)
_nii2zarr_prepare_out(x::Zarr.ZarrGroup) = x
_nii2zarr_prepare_chunk(x::NTuple{N, Integer} where N) = (x[begin:3],)
_nii2zarr_prepare_chunk(x::Integer) = ((x, x, x, 1, 1),)
_nii2zarr_prepare_chunk(x::NTuple{1, Integer}) = ((x..., x..., x..., 1, 1),)
_nii2zarr_prepare_chunk(x::NTuple{2, Integer}) = ((x..., x[end], 1, 1),)
_nii2zarr_prepare_chunk(x::NTuple{3, Integer}) = ((x, 1, 1),)
_nii2zarr_prepare_chunk(x::NTuple{4, Integer}) = ((x, 1),)
_nii2zarr_prepare_chunk(x::NTuple{5, Integer}) = (x,)
_nii2zarr_prepare_chunk(x::Tuple{Tuple}) = _nii2zarr_prepare_chunk(collect(x))
_nii2zarr_prepare_chunk(x::Array) = map(_nii2zarr_prepare_chunk, x)
_nii2zarr_prepare_chunk(x::Array{NTuple{5, Integer}}) = x
_nii2zarr_prepare_label(x) = Zarr.zopen(x)
_nii2zarr_prepare_label(x::Zarr.ZarrGroup) = x
_nii2zarr_prepare_compressor(x::String) = getfield(CompressorType, Symbol(x))
_nii2zarr_prepare_compressor(x::Symbol) = getfield(CompressorType, x)
_nii2zarr_prepare_compressor(x::CompressorType) = x
_nii2zarr_prepare_method(x::String) = getfield(MethodType, Symbol(x))
_nii2zarr_prepare_method(x::Symbol) = getfield(MethodType, x)
_nii2zarr_prepare_method(x::MethodType) = x


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

_ensure5d(::AbstractArray{T, N} where N) = error("Too many dimensions for conversion to nii.zarr")
_ensure5d(x::AbstractArray{T, 5} where T) = x
_ensure5d(x::AbstractArray{T, 4} where T) = reshape(x, (size(x)..., 1))
_ensure5d(x::AbstractArray{T, 3} where T) = reshape(x, (size(x)..., 1, 1))
_ensure5d(x::AbstractArray{T, 2} where T) = reshape(x, (size(x)..., 1, 1, 1))
_ensure5d(x::AbstractArray{T, 1} where T) = reshape(x, (size(x)..., 1, 1, 1, 1))
_ensure5d(x::AbstractArray{T, 0} where T) = reshape(x, (1, 1, 1, 1, 1))


function _nii2zarr_ome_attrs(shapes, niiheader)
    attrs = Dict(
        "nifti" => niiheader,
        "multiscales" => [Dict()]
    )
    multiscales = attrs["multiscales"][1]

    multiscales["axes"] = [
        Dict(
            "name" => "t",
            "type" => "time",
            "unit" => niiheader["units"]["time"],
        ),
        Dict(
            "name" => "c",
            "type" => "channel"
        ),
        Dict(
            "name" => "z",
            "type" => "space",
            "unit" => niiheader["units"]["space"],
        ),
        Dict(
            "name" => "y",
            "type" => "space",
            "unit" => niiheader["units"]["space"],
        ),
        Dict(
            "name" => "x",
            "type" => "space",
            "unit" => niiheader["units"]["space"],
        )
    ]

    multiscales["datasets"] = empty(Vector{Dict}, length(shapes))
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
                "scale" => [
                    1.0,
                    1.0,
                    (shapes[0][0]/shapes[n][0])*jsonheader["pixdim"][2],
                    (shapes[0][1]/shapes[n][1])*jsonheader["pixdim"][1],
                    (shapes[0][2]/shapes[n][2])*jsonheader["pixdim"][0],
                ]
            ),
            Dict(
                "type" => "translation",
                "translation" => [
                    1.0,
                    1.0,
                    (shapes[0][0]/shapes[n][0] - 1)*jsonheader["pixdim"][2]*0.5,
                    (shapes[0][1]/shapes[n][1] - 1)*jsonheader["pixdim"][1]*0.5,
                    (shapes[0][2]/shapes[n][2] - 1)*jsonheader["pixdim"][0]*0.5,
                ]
            ),
        ]
    end

    # time scale
    multiscales["coordinateTransformations"] = [Dict(
        "scale" => [
            len(niiheader["pixdim"]) > 3 ? niiheader["pixdim"][3] : 1.0,
            1.0,
            1.0,
            1.0,
            1.0
        ],
        "type" => "scale"
    )]

    return attrs
end


nii2zarr(inp.NIVolume, out;
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
  `NIfTI.read` handles.
- The output should be a `Zarr.ZarrGroup`, but can also be anything that
  `Zarr.zgroup` or `Zarr.zopen` handle.

# Keyword arguments
- `chunk`: Chunk size, per x/y/z/t/c, per level.
    * The inner tuple allows different chunk sizes to be used along
      each dimension.
    * The outer array allows different chunk sizes to be used at
      different pyramid levels.
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
function nii2zarr(inp::NIfTI.NIVolume, out::Zarr.ZarrGroup;
    chunk::Array{NTuple{5, Integer}} = [(64, 64, 64, 1, 1)],
    nb_levels::Integer = -1,
    method::PyramidMethod = PyramidMethod.gaussian,
    label::IsLabel = IsLabel.auto,
    fill_value::Real = nothing,
    compressor::CompressorType = CompressorType.blosc,
    compressor_options...
)
    # Convert NIfTI.jl's header into my header (they only handles nifti-1)
    header = reinterpret(Nifti1Header, header(inp))
    jsonheader = nii2json(header)

    # Fix array shape
    data = _ensure5d(inp.raw)
    data = PermutedDimsArray(data, [4, 5, 3, 2, 1])

    # Compute image pyramid
    bool_label = jsonheader["intent"]["code"] in ("LABEL", "NEURONAMES")
    if label == IsLabel.yes
        bool_label = true
    elseif label == IsLabel.no
        bool_label = false
    end

    data = _make_pyramid5d(data, nb_levels, label)
    nb_levels = len(data)
    shapes = [d.shape[3:end] for d in data]
    attrs = _nii2zarr_ome_attrs(shapes, jsonheader)

    # Prepare array metadata at each level
    compressor = _make_compressor(compressor; compressor_options...)
    chunk = _nii2zarr_prepare_chunk(chunk, nb_levels)
    chunk = map(x->(x[4], x[5], x[3], x[2], x[1]), chunk)
    chunk = [Dict(
        "chunks" => c,
        # "dimension_separator" => '/',
        # "order" => 'F',
        # "dtype" => jsonheader["datatype"],
        "fill_value" => fill_value,
        "compressor" => compressor,
    ) for c in chunk]

    # Write group attributes
    writeattrs(out.storage, out.path, attrs)

    # Write zarr arrays
    for n in eachindex(shapes)
        subarray = Zarr.zcreate(
            eltype(data[n]), out, string(n-1), shapes[n]...;
            attrs = Dict("_ARRAY_DIMENSIONS" => ["time", "channel", "z", "y", "x"])
            chunk[n]...,
        )
        subarray[:,:,:,:,:] = data[n]
    end

end
