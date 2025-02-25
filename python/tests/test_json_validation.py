import json
import os.path as op
import tempfile
import unittest

import nibabel as nib
from jsonschema import validate
from jsonschema.exceptions import ValidationError

import niizarr

HERE = op.dirname(op.abspath(__file__))
DATA = op.join(HERE, "data")
ROOT = op.join(HERE, "..", "..")


class TestJSONValidation(unittest.TestCase):

    def setUp(self):
        self.schema_file = op.join(ROOT, "nifti-zarr-schema-0.3.json")
        with open(self.schema_file) as f:
            self.schema = json.load(f)

    def test_json_validation(self):
        test_files = ["example_nifti2.nii.gz", "example4d.nii.gz"]
        for nifti_file in test_files:
            with self.subTest(nifti_file=nifti_file):
                data = nib.load(op.join(DATA, nifti_file))
                with tempfile.TemporaryDirectory() as tmpdir:
                    zarr_file = op.join(tmpdir, "test.nii.zarr")
                    niizarr.nii2zarr(data, zarr_file)

                    json_file = op.join(zarr_file, "nifti/.zattrs")
                    with open(json_file) as f:
                        json_obj = json.load(f)

                    try:
                        validate(json_obj, self.schema)
                    except ValidationError:
                        self.fail(json.dumps(json_obj))


if __name__ == '__main__':
    unittest.main()
