import gzip
import tempfile
import unittest

import nibabel as nib
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
