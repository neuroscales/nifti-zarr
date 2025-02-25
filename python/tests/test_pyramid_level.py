import os
import tempfile
import unittest

import numpy as np
from nibabel import Nifti1Image

from niizarr import nii2zarr


def count_levels(nifti_zarr):
    """
    In a .nii.zarr path, the number of folders should be the number of the levels plus 1(nifti folder)
    """
    for root, dir, file in os.walk(nifti_zarr):
        return len(dir) - 1


class TestPyramidLevel(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_zarr = os.path.join(self.temp_dir.name, 'output.nii.zarr')
        img = np.random.rand(16, 32, 64)
        self.ni = Nifti1Image(img, np.eye(4))

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_no_nb_level(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=1)
        self.assertEqual(1, count_levels(self.output_zarr))

    def test_constant_nb_level(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=2)
        self.assertEqual(2, count_levels(self.output_zarr))

    def test_chunk_size_16(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=-1, chunk=16)
        self.assertEqual(3, count_levels(self.output_zarr))

    def test_chunk_size_24(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=-1, chunk=24)
        self.assertEqual(3, count_levels(self.output_zarr))

    def test_chunk_size_32(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=-1, chunk=32)
        self.assertEqual(2, count_levels(self.output_zarr))

    def test_chunk_size_48(self):
        # 64 -> 32 OK
        nii2zarr(self.ni, self.output_zarr, nb_levels=-1, chunk=48)
        self.assertEqual(2, count_levels(self.output_zarr))

    def test_chunk_size_64(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=-1, chunk=64)
        self.assertEqual(1, count_levels(self.output_zarr))

    def test_chunk_size_96(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=-1, chunk=96)
        self.assertEqual(1, count_levels(self.output_zarr))

    def test_chunk_size_512(self):
        nii2zarr(self.ni, self.output_zarr, nb_levels=-1, chunk=512)
        self.assertEqual(1, count_levels(self.output_zarr))


if __name__ == '__main__':
    unittest.main()
