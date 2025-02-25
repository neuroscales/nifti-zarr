import os.path as op
import tempfile
import unittest

import nibabel as nib
import numpy as np

from niizarr import zarr2nii

HERE = op.dirname(op.abspath(__file__))
DATA = op.join(HERE, "data")


class Testzarr2nii(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_same_result_nifti1(self):
        zarr_file = op.join(DATA, "example4d.nii.zarr")
        nifti_file = op.join(DATA, "example4d.nii.gz")

        converted = zarr2nii(zarr_file)
        loaded = nib.load(nifti_file)

        self.assertEqual(str(loaded.header), str(converted.header))
        self.assertEqual(str(loaded.header.extensions),
                         str(converted.header.extensions))
        np.testing.assert_array_almost_equal(loaded.get_fdata(), converted.get_fdata())

    def test_same_result_nifti2(self):
        zarr_file = op.join(DATA, "example_nifti2.nii.zarr")
        nifti_file = op.join(DATA, "example_nifti2.nii.gz")

        converted = zarr2nii(zarr_file)
        loaded = nib.load(nifti_file)

        self.assertEqual(str(loaded.header), str(converted.header))
        self.assertEqual(str(loaded.header.extensions),
                         str(converted.header.extensions))
        np.testing.assert_array_almost_equal(loaded.get_fdata(), converted.get_fdata())
