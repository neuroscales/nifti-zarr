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
    except:
        print("faile to load sample data from skimage, using random data")
        img = np.random.rand(100, 200, 300)
    ni = klass_map[version](img, np.eye(4))

    return ni
