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

nii2json(header::NIfTI.NIfTIHeader) = nii2json(header, false)

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
function nii2json(header::NIfTI.NIfTIHeader, has_extensions::Bool)

    magic = _tuple2string(header.magic)
 
    magic = magic[1] * 'i' * magic[3]
    ndim = header.dim[1]
    
    intent_code = IntentRecoder(header.intent_code)

    nb_intent_prm = IntentNbPrm[IntentRecoder(header.intent_code)]

    jsonheader = Dict(
        "NIIFormat" => magic,
        "Dim" => collect(header.dim[2:2+ndim-1]),
        "VoxelSize" => collect(header.pixdim[2:2+ndim-1]),
        "Unit" => Dict(
            "L" => UnitRecoder(header.xyzt_units & 0x07),
            "T" => UnitRecoder(header.xyzt_units & 0x38)
        ),
        "DataType" => DataTypeRecoder(header.datatype),
        "DimInfo" => Dict(
            "Freq" => header.dim_info & 0x03,
            "Phase" => (header.dim_info >> 2) & 0x03,
            "Slice" => (header.dim_info >> 4) & 0x03
        ),
        "Intent" => intent_code,
        "Name" => _tuple2string(header.intent_name),
        "Param1" => nothing,
        "Param2" => nothing,
        "Param3" => nothing,
        "ScaleSlope"=> header.scl_slope,
        "ScaleOffset"=> header.scl_inter,
        "FirstSliceID" => header.slice_start,
        "LastSliceID"=> header.slice_end,
        "SliceType" =>  SliceOrderRecoder(header.slice_code),
        "SliceTime" => header.slice_duration,
        "MinIntensity" => header.cal_min,
        "MaxIntensity" => header.cal_max,
        "TimeOffset" => header.toffset,
        "Description" => _tuple2string(header.descrip),
        "AuxFile" => _tuple2string(header.aux_file),
        "QForm"=> XFormRecoder(header.qform_code),
        "Quatern"=> Dict(
            "b"=> header.quatern_b,
            "c"=> header.quatern_c,
            "d"=> header.quatern_d,
        ),
        "QuaternOffset"=> Dict(
            "x"=> header.qoffset_x,
            "y"=> header.qoffset_y,
            "z"=> header.qoffset_z,
        ),
        "Orientation"=> Dict(
            "x"=> header.pixdim[1] == 0 ? "r" : "l",
            "y"=> "a",
            "z"=> "s",
        ),
        "SForm"=> XFormRecoder(header.sform_code),
        "Affine" => [[header.srow_x...] [header.srow_y...] [header.srow_z...]],
        "NIFTIExtension" => [has_extensions ? 1 : 0, 0 , 0 , 0]
    )
    if !isfinite(jsonheader["ScaleSlope"])
        jsonheader["ScaleSlope"] = 0.0
    end
    if !isfinite(jsonheader["ScaleOffset"])
        jsonheader["ScaleOffset"] = 0.0
    end

    for (i,v) in pairs(nb_intent_prm)
        jsonheader["Param$i"] = v
    end


    return jsonheader
end


function _make_pyramid3d(data3d::AbstractArray, nb_levels::Integer, label::Bool=false)
    nbatch = size(data3d)[1:end-3]
    nxyz = size(data3d)[end-2:end] 
    data4d = reshape(data3d, (prod(nbatch), nxyz...))
    nb_levels = max(nb_levels, 1)
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
    
    if eltype(data3d) != eltype(levels[1])
        if eltype(data3d)<: Integer && eltype(levels[1])<: AbstractFloat
            levels = map(x->round.(x), levels)
        end
        levels = map(x->convert(Array{eltype(data3d)},x), levels)
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
    # TODO: missing attrs that exists in python, hardcoded now
    multiscales["version"] = "0.4"
    multiscales["name"] = "/"
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
        axis == "c"
        ? Dict(
            "name" => axis,
            "type" => typemap[axis]
        )
        : Dict(
            "name" => axis,
            "type" => typemap[axis],
            "unit" => typemap[axis]=="time" ? JNIfTIUnit2ZarrUnit[niiheader["Unit"]["T"]] : JNIfTIUnit2ZarrUnit[niiheader["Unit"]["L"]]
        )
    ) for axis in axes]

    multiscales["datasets"] = [Dict() for _ in 1:length(shapes)]
    for n in eachindex(shapes)
        # spatial scale at each pyramid level
        level = multiscales["datasets"][n]

        level["path"] = string(n-1)
        # Image.gaussian_pyramid aligns voxel edges across scales
        # (same behavior as scipy.ndimage.zoom(..., grid_mode=True))
        # so the effective scaling is the shape ratio, and there is
        # a half voxel shift wrt to the "center of first voxel" frame
        level["coordinateTransformations"] = [
            Dict(
                "type" => "scale",
                "scale" => [(
                    typemap[axis] == "space"
                    ? (shapes[1][i]/shapes[n][i])*niiheader["VoxelSize"][idxmap[axis]]
                    : 1.0
                ) for (i, axis) in enumerate(axes)]
            ),
            Dict(
                "type" => "translation",
                "translation" => [(
                    typemap[axis] == "space"
                    ? (shapes[1][i]/shapes[n][i] - 1)*niiheader["VoxelSize"][idxmap[axis]]*0.5
                    : 0.0
                ) for (i, axis) in enumerate(axes)]
            ),
        ]
    end

    # time scale
    multiscales["coordinateTransformations"] = [Dict(
        "scale" => [(
            typemap[axis] == "time"
            ? niiheader["VoxelSize"][idxmap[axis]]
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
    chunk::Vector{NTuple{5, T}} where T <: Integer = [(128, 128, 128, 1, 1)],
    nb_levels::Integer = -1,
    method::PyramidMethod.T = PyramidMethod.gaussian,
    label::IsLabel.T = IsLabel.auto,
    fill_value = nothing,
    compressor::CompressorType.T = CompressorType.blosc,
    compressor_options...
)
    header = inp.header
    jsonheader = nii2json(header, !isempty(inp.extensions))
    
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
    
    # set bool label based on override parameter or from header
    bool_label = label == IsLabel.yes || (label == IsLabel.auto && jsonheader["Intent"] in ("label", "neuronames"))

    # Compute image pyramid
    if nb_levels == -1
        nxyz = collect(size(data)[end-2:end])
        spacial_chunksize = collect(chunk[1][1:3])
        default_nb_levels = ceil(Int, log2(maximum(nxyz ./ spacial_chunksize))) + 1
        # make sure at lease 1 level
        default_nb_levels = max(default_nb_levels, 1)
        nb_levels=default_nb_levels
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
    binheader = bytesencode(header)
    # NIfTI.jl has problem with handling extensions, so we cannot write it for now.
    # if length(inp.extensions) > 0
    #     io = IOBuffer()
    #     write(io, inp.extensions)
    #     seekstart(io)          # Reset the position of the buffer to read from the beginning
    #     ext_bytes = read(io)   # Read the buffer content into a UInt8 vector
    #     binheader = vcat(binheader, ext_bytes)
    # end

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
    
    # Remove all exisitng levels
    for k in keys(out.arrays)
        if isa(tryparse(Float64,k), Number)
            delete!(out.storage, "", k)
        end
    end

    # Write zarr arrays
    for n in eachindex(shapes)
        # everything is reversed since zarr.jl write things in reversed order, 
        # see https://github.com/JuliaIO/Zarr.jl/issues/78
        subarray = Zarr.zcreate(
            eltype(data[n]), out, string(n-1), reverse(shapes[n])...;
            chunks=reverse(chunk[n]),
            fill_as_missing=false,
            fill_value=nothing,
            compressor=compressor,
            # TODO: dimension separator and order is not yet supported in Zarr.jl
            # dimension_separator='/',
            # order='F',
            attrs = Dict("_ARRAY_DIMENSIONS" => ARRAY_DIMENSIONS)
        )
        # subarray[:,:,:] = data[n]
        copy!(subarray, PermutedDimsArray(data[n], reverse(1:ndims(data[n]))))
    end
    
end

# # TODO: Move this to nii2zarr
#     # Fix data type
#     bo = ByteOrderSymbol[byteorder]
#     if isa(jsonheader["DataType"], Array)
#         jsonheader["DataType"] = map(
#             x -> [x[1], '|' * x[2]],
#             jsonheader["DataType"]
#         )
#     elseif jsonheader["DataType"][end] == '1'
#         jsonheader["DataType"] = '|' + jsonheader["DataType"]
#     else
#         jsonheader["DataType"] = bo * jsonheader["DataType"]
#     end

# Fix delete! from Zarr.jl, which implements it as:
#   Base.delete!(s::DirectoryStore, k::String) = isfile(joinpath(s.folder, k)) && rm(joinpath(s.folder, k))
# Their version only deletes files, and not directories.
delete!(s::Zarr.DirectoryStore, k::String) = rm(joinpath(s.folder, k); force=true, recursive=true)
