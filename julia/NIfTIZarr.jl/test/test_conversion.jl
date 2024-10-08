using Test
using NIfTIZarr
import NIfTI
using JSON
using JSONSchema
import Base: ==, hash

function ==(a::NIfTI.NIfTI1Header, b::NIfTI.NIfTI1Header)
    # Check if both objects are of the same type
    if typeof(a) !== typeof(b)
        return false
    end
    # Iterate over all field names
    for name in fieldnames(NIfTI.NIfTI1Header)
        field_a = getfield(a, name)
        field_b = getfield(b, name)
        if field_a != field_b
            return false
        end
    end
    return true
end


original_file = joinpath(dirname(@__FILE__), "data/example4d.nii.gz")
const TEMP_DIR_NAME = mktempdir()
output_file = joinpath(TEMP_DIR_NAME, "$(tempname()).nii.zarr")

NIfTIZarr.nii2zarr(original_file, output_file, chunk=128)
back = NIfTIZarr.zarr2nii(output_file)
original = NIfTI.niread(original_file)


@test isequal(original.header, back.header)
schema = Schema(JSON.parsefile(joinpath(dirname(@__FILE__), "data/nifti-zarr-schema-0.3.json")))
output_json =  JSON.parsefile(joinpath(output_file, "nifti/.zattrs"))
# println(output_json["QForm"])
@test isvalid(schema, output_json)
# validate(schema, output_json)
# @test isvalid()
