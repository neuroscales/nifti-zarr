import filecmp
import json
import os

import numpy as np
import skimage as sk
from nibabel import Nifti1Image, Nifti2Image

klass_map = {
    1: Nifti1Image,
    2: Nifti2Image
}


def get_nifti_image(version=1):
    try:
        img = sk.data.brain()
    except Exception as e:
        print(e)
        print("failed to load sample data from skimage, using random data")
        img = np.random.rand(100, 200, 300)
    ni = klass_map[version](img, np.eye(4))

    return ni


class dircmp(filecmp.dircmp):
    """
    Compare the content of dir1 and dir2. In contrast with filecmp.dircmp, this
    subclass compares the content of files with the same path.
    """

    def phase3(self):
        """
        Find out differences between common files.
        Ensure we are using content comparison with shallow=False.
        """
        fcomp = filecmp.cmpfiles(self.left, self.right, self.common_files,
                                 shallow=True)
        self.same_files, self.diff_files, self.funny_files = fcomp


def is_same(dir1, dir2):
    """
    Compare two directory trees content.
    Return False if they differ, True is they are the same.
    """
    compared = dircmp(dir1, dir2)
    if (compared.left_only or compared.right_only or compared.diff_files
            or compared.funny_files):
        return False
    for subdir in compared.common_dirs:
        if not is_same(os.path.join(dir1, subdir), os.path.join(dir2, subdir)):
            return False
    return True


def compare_json_objects(obj1, obj2):
    """Recursively compare two JSON objects."""
    if isinstance(obj1, dict) and isinstance(obj2, dict):
        if obj1.keys() != obj2.keys():
            return False
        for key in obj1:
            if not compare_json_objects(obj1[key], obj2[key]):
                return False
    elif isinstance(obj1, list) and isinstance(obj2, list):
        if len(obj1) != len(obj2):
            return False
        for i in range(len(obj1)):
            if not compare_json_objects(obj1[i], obj2[i]):
                return False
    else:
        if obj1 != obj2:
            return False
    return True


def compare_zarr_archives(path1, path2):
    """
    Compare two Zarr archives by comparing the `.zarray`, `.zattrs`, and `nifti/0` JSON files
    in the root and subdirectories.

    :param path1: Path to the first Zarr archive.
    :param path2: Path to the second Zarr archive.
    :return: True if both archives are identical, False otherwise.
    """
    # Traverse both directories
    for root, _, files in os.walk(path1):
        relative_path = os.path.relpath(root, path1)
        other_root = os.path.join(path2, relative_path)

        # Compare relevant files
        for file_name in ['.zarray', '.zattrs']:
            file_path1 = os.path.join(root, file_name)
            file_path2 = os.path.join(other_root, file_name)
            # Check if both files exist
            if os.path.exists(file_path1) and os.path.exists(file_path2):
                with open(file_path1, 'r') as f1, open(file_path2, 'r') as f2:
                    try:
                        json1 = json.load(f1)
                        json2 = json.load(f2)
                    except json.JSONDecodeError as e:
                        print(f"Error decoding JSON in {file_path1} or {file_path2}: {e}")
                        return False

                    if not compare_json_objects(json1, json2):
                        print(f"Mismatch found in {file_name} at {relative_path}")
                        return False
            else:
                # If one file exists and the other does not
                if os.path.exists(file_path1) or os.path.exists(file_path2):
                    print(f"File missing in one of the directories: {file_path1} or {file_path2}")
                    return False

    return True

# This script generates trusted test data.
if __name__ == '__main__':
    from niizarr import *
    import zarr

    input_files = ["data/example4d.nii.gz", "data/example_nifti2.nii.gz"]
    for input_file in input_files:
        output_file = input_file.replace(".gz", ".zarr")
        nii2zarr(input_file, output_file, chunk=64)
        inp = zarr.storage.DirectoryStore(output_file)

        inp = zarr.group(store=inp)
        for layer in (0, 1, 'nifti'):
            if str(layer) not in inp:
                continue
            content = np.array(inp[str(layer)])

            np.save(os.path.join(output_file, f"{layer}.npy"), content)
