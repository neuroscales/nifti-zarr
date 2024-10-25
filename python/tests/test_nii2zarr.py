import gzip
import os
import tempfile
import unittest

import nibabel as nib
import numpy as np
import zarr
from _data import compare_zarr_archives
from niizarr import nii2zarr


class Testnii2zarr(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_non_exist_file(self):
        self.assertRaises(Exception, nii2zarr, 'non_exist_file.nii', self.temp_dir.name)

    def test_input_path(self):
        nii2zarr("data/example4d.nii.gz", self.temp_dir.name)

    def test_input_fd(self):
        with gzip.open("data/example4d.nii.gz", "rb") as f:
            nii2zarr(f, self.temp_dir.name)

    def test_input_Nifti1Image(self):
        ni = nib.load("data/example4d.nii.gz")
        nii2zarr(ni, self.temp_dir.name)

    def test_input_Nifti2Image(self):
        ni = nib.load("data/example_nifti2.nii.gz")
        nii2zarr(ni, self.temp_dir.name)

    def test_same_result_nifti1(self):
        written_zarr = os.path.join(self.temp_dir.name, "example4d.nii.zarr")
        nii2zarr("data/example4d.nii.gz", written_zarr, chunk=64)
        self.assertTrue(compare_zarr_archives(written_zarr, "data/example4d.nii.zarr"))
        written_data = zarr.load(written_zarr)
        for layer in (0, 1, 'nifti'):
            reference_data = np.load(f"data/example4d.nii.zarr/{layer}.npy")
            np.testing.assert_array_almost_equal(reference_data, written_data[layer])

    def test_same_result_nifti2(self):
        written_zarr = os.path.join(self.temp_dir.name, "example_nifti2.nii.zarr")
        nii2zarr("data/example_nifti2.nii.gz", written_zarr, chunk=64)
        self.assertTrue(compare_zarr_archives(written_zarr, "data/example_nifti2.nii.zarr"))
        written_data = zarr.load(written_zarr)
        for layer in (0, 'nifti'):
            reference_data = np.load(f"data/example_nifti2.nii.zarr/{layer}.npy")
            np.testing.assert_array_almost_equal(reference_data, written_data[layer])
