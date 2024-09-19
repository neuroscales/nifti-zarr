import unittest
import nibabel as nib
import tempfile
import os
import json
from jsonschema import validate
from jsonschema.exceptions import ValidationError

from python import niizarr


class TestNiizarrConversion(unittest.TestCase):

    def setUp(self):
        self.schema_file = "../../nifti-zarr-schema-0.3.json"
        with open(self.schema_file) as f:
            self.schema = json.load(f)

    def test_json_validation(self):
        test_files = ["data/example_nifti2.nii.gz", "data/example4d.nii.gz"]
        for nifti_file in test_files:
            with self.subTest(nifti_file=nifti_file):
                data = nib.load(nifti_file)
                with tempfile.TemporaryDirectory() as tmpdir:
                    zarr_file = os.path.join(tmpdir, "test.nii.zarr")
                    niizarr.nii2zarr(data, zarr_file)

                    json_file = os.path.join(zarr_file, "nifti/.zattrs")
                    with open(json_file) as f:
                        json_obj = json.load(f)

                    try:
                        validate(json_obj, self.schema)
                    except ValidationError as e:
                        self.fail(json.dumps(json_obj))




if __name__ == '__main__':
    unittest.main()
