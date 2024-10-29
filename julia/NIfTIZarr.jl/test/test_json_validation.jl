using Test
using NIfTIZarr
import NIfTI
using JSON
using JSONSchema

schema = Schema(JSON.parsefile(joinpath(dirname(@__FILE__), "data/nifti-zarr-schema-0.3.json")));

function test_json_validation(input_file::String)
    original_file = input_file
    output_file = "$(tempname()).nii.zarr";
    NIfTIZarr.nii2zarr(original_file, output_file, chunk=128);
    output_json = JSON.parsefile(joinpath(output_file, "nifti/.zattrs"));

    @test isvalid(schema, output_json)
end


test_json_validation(joinpath(dirname(@__FILE__), "data/example4d.nii.gz"))
test_json_validation(joinpath(dirname(@__FILE__), "data/example_nifti2.nii.gz"))

