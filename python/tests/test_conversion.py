import unittest
import nibabel as nib
import tempfile
import os

from python import niizarr


class TestNiizarrConversion(unittest.TestCase):

    def test_conversion_roundtrip_header(self):
        test_files = ["data/example_nifti2.nii.gz", "data/example4d.nii.gz"]
        for nifti_file in test_files:
            with self.subTest(nifti_file=nifti_file):
                data = nib.load(nifti_file)
                with tempfile.TemporaryDirectory() as tmpdir:
                    zarr_file = os.path.join(tmpdir, "test.ome.zarr")

                    niizarr.nii2zarr(data, zarr_file)

                    loaded = niizarr.zarr2nii(zarr_file)

                    original_header = data.header
                    loaded_header = loaded.header

                    self.assertEqual(str(original_header), str(loaded_header))

    def test_conversion_roundtrip_extension(self):
        nifti_file = "data/example_nifti2.nii.gz"
        data = nib.load(nifti_file)
        with tempfile.TemporaryDirectory() as tmpdir:
            zarr_file = os.path.join(tmpdir, "test.ome.zarr")

            niizarr.nii2zarr(data, zarr_file)

            loaded = niizarr.zarr2nii(zarr_file)

            original_header = data.header
            loaded_header = loaded.header

            self.assertEqual(str(original_header.extensions), str(loaded_header.extensions))
if __name__ == '__main__':
    unittest.main()


"""
sizeof_hdr      : 348
data_type       : b''
db_name         : b''
extents         : 0
session_error   : 0
regular         : b'r'
dim_info        : 57
dim             : [  4 128  96  24   2   1   1   1]
intent_p1       : 0.0
intent_p2       : 0.0
intent_p3       : 0.0
intent_code     : none
datatype        : int16
bitpix          : 16
slice_start     : 0
pixdim          : [-1.000000e+00  2.000000e+00  2.000000e+00  2.199999e+00  2.000000e+03
  1.000000e+00  1.000000e+00  1.000000e+00]
vox_offset      : 0.0
scl_slope       : nan
scl_inter       : nan
slice_end       : 23
slice_code      : unknown
xyzt_units      : 10
cal_max         : 1162.0
cal_min         : 0.0
slice_duration  : 0.0
toffset         : 0.0
glmax           : 0
glmin           : 0
descrip         : b'FSL3.3\x00 v2.25 NIfTI-1 Single file format'
aux_file        : b''
qform_code      : scanner
sform_code      : scanner
quatern_b       : -1.9451068e-26
quatern_c       : -0.9967085
quatern_d       : -0.08106874
qoffset_x       : 117.8551
qoffset_y       : -35.722942
qoffset_z       : -7.2487984
srow_x          : [-2.0000000e+00  6.7147157e-19  9.0810245e-18  1.1785510e+02]
srow_y          : [-6.7147157e-19  1.9737115e+00 -3.5552824e-01 -3.5722942e+01]
srow_z          : [ 8.2554809e-18  3.2320762e-01  2.1710818e+00 -7.2487984e+00]
intent_name     : b''
magic           : b'n+1'

"""