import io
import json
import os
import tempfile
import unittest

import nibabel as nib
import numpy as np
from nibabel import Nifti1Header, Nifti1Image
from nibabel.nifti1 import Nifti1Extension

from _data import get_nifti_image
from niizarr import nii2zarr
from niizarr._header import SYS_BYTEORDER, HEADERTYPE1


class TestEndian(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.sys_endian_nifti = os.path.join(self.tempdir.name, "sys_endian.nii")
        self.inv_endian_nifti = os.path.join(self.tempdir.name, "inv_endian.nii")
        self.sys_endian_zarr = os.path.join(self.tempdir.name, "sys_endian.nii.zarr")
        self.inv_endian_zarr = os.path.join(self.tempdir.name, "inv_endian.nii.zarr")

        ni = get_nifti_image()
        # so there will be a endian sign in dtype
        ni.set_data_dtype(np.dtype("int32"))

        header = np.frombuffer(ni.header.structarr.tobytes(), HEADERTYPE1).byteswap().newbyteorder()
        header_file_obj = io.BytesIO(header.tobytes())
        header = Nifti1Header.from_fileobj(header_file_obj)
        inv_ni = Nifti1Image(ni.get_fdata(), np.eye(4), header=header)

        extension = Nifti1Extension(0, b"This is some extension content")
        ni.header.extensions.append(extension)
        inv_ni.header.extensions.append(extension)

        nib.save(ni, self.sys_endian_nifti)
        nib.save(inv_ni, self.inv_endian_nifti)

        nii2zarr(self.sys_endian_nifti, self.sys_endian_zarr)
        nii2zarr(self.inv_endian_nifti, self.inv_endian_zarr)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_setup_endian(self):
        self.assertTrue(nib.load(self.sys_endian_nifti).header.endianness == SYS_BYTEORDER)
        self.assertFalse(nib.load(self.inv_endian_nifti).header.endianness == SYS_BYTEORDER)

    def test_system_endian_header(self):
        with open(self.sys_endian_nifti, "rb") as f:
            sys_nifti_content = f.read(4)
        with open(os.path.join(self.sys_endian_zarr, "nifti/0"), "rb") as f:
            sys_zarr_content = f.read(4)
        self.assertEqual(sys_nifti_content, sys_nifti_content)

    def test_inv_endian_header(self):
        with open(self.inv_endian_nifti, "rb") as f:
            inv_nifti_content = f.read(4)
        with open(os.path.join(self.inv_endian_zarr, "nifti/0"), "rb") as f:
            inv_zarr_content = f.read(4)
        self.assertEqual(inv_nifti_content, inv_zarr_content)

    def test_system_endian_zarr_content(self):
        with open(os.path.join(self.sys_endian_zarr, "0/.zarray"), "r") as f:
            zarray = json.load(f)
        self.assertEqual(zarray["dtype"][0], SYS_BYTEORDER)

    def test_inv_endian_zarr_content(self):
        with open(os.path.join(self.inv_endian_zarr, "0/.zarray"), "r") as f:
            zarray = json.load(f)
        self.assertNotEqual(zarray["dtype"][0], SYS_BYTEORDER)

    def test_system_endian_extension(self):

        with open(os.path.join(self.sys_endian_zarr, "nifti/0"), "rb") as f:
            f.seek(348)
            zarr_extension_content = f.read(4)
        with open(self.sys_endian_nifti, "rb") as f:
            f.seek(348 + 4)
            nifti_extension_content = f.read(4)
        self.assertEqual(zarr_extension_content, nifti_extension_content)
        if SYS_BYTEORDER == "<":
            self.assertEqual(zarr_extension_content, b'\x30\x00\x00\x00')
        else:
            self.assertEqual(zarr_extension_content, b'\x00\x00\x00\x30')

    def test_inv_endian_extension(self):

        with open(os.path.join(self.inv_endian_zarr, "nifti/0"), "rb") as f:
            f.seek(348)
            zarr_extension_content = f.read(4)
        with open(self.inv_endian_nifti, "rb") as f:
            f.seek(348 + 4)
            nifti_extension_content = f.read(4)
        self.assertEqual(zarr_extension_content, nifti_extension_content)
        if SYS_BYTEORDER == ">":
            self.assertEqual(zarr_extension_content, b'\x30\x00\x00\x00')
        else:
            self.assertEqual(zarr_extension_content, b'\x00\x00\x00\x30')
