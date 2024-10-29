using Test
using Random
using LinearAlgebra
import NIfTI
using Zarr
using NIfTIZarr

function count_levels(nifti_zarr::String)
    """
    In a .nii.zarr path, the number of folders should be the number of the levels plus 1(nifti folder)
    """
    for (root, dirs, files) in walkdir(nifti_zarr)
        return length(dirs) - 1
    end
end

temp_dir = mktempdir()
output_zarr = joinpath(temp_dir, "output.nii.zarr")
img_data = rand(Float64,(16, 32, 64))
ni = NIfTI.NIVolume(img_data)

@testset "No NB Level" begin
    nii2zarr(ni, output_zarr, nb_levels=1)
    @test count_levels(output_zarr) == 1
end

@testset "Constant NB Level" begin
    nii2zarr(ni, output_zarr, nb_levels=2)
    @test count_levels(output_zarr) == 2
end

@testset "Chunk Size 16" begin
    nii2zarr(ni, output_zarr, nb_levels=-1, chunk=16)
    @test count_levels(output_zarr) == 3
end

@testset "Chunk Size 24" begin
    nii2zarr(ni, output_zarr, nb_levels=-1, chunk=24)
    @test count_levels(output_zarr) == 3
end

@testset "Chunk Size 32" begin
    nii2zarr(ni, output_zarr, nb_levels=-1, chunk=32)
    @test count_levels(output_zarr) == 2
end

@testset "Chunk Size 48" begin
    nii2zarr(ni, output_zarr, nb_levels=-1, chunk=48)
    @test count_levels(output_zarr) == 2
end

@testset "Chunk Size 64" begin
    nii2zarr(ni, output_zarr, nb_levels=-1, chunk=64)
    @test count_levels(output_zarr) == 1
end

@testset "Chunk Size 96" begin
    nii2zarr(ni, output_zarr, nb_levels=-1, chunk=96)
    @test count_levels(output_zarr) == 1
end

@testset "Chunk Size 512" begin
    nii2zarr(ni, output_zarr, nb_levels=-1, chunk=512)
    @test count_levels(output_zarr) == 1
end
