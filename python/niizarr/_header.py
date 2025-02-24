import re
import sys
import warnings

import numpy as np
from nibabel import (Nifti1Image, Nifti1Header, Nifti2Image, Nifti2Header)

SYS_BYTEORDER = '<' if sys.byteorder == 'little' else '>'
SYS_BYTEORDER_SWAPPED = '>' if SYS_BYTEORDER == '<' else '<'

NIFTI_1_HEADER_SIZE = 348
NIFTI_2_HEADER_SIZE = 540

HEADERTYPE1 = np.dtype([
    ("sizeof_hdr", "i4"),
    ("data_type", "S10"),
    ("db_name", "S18"),
    ("extents", "i4"),
    ("session_error", "i2"),
    ("regular", "u1"),
    ("dim_info", "u1"),
    ("dim", "i2", 8),
    ("intent_p", "f4", 3),
    ("intent_code", "i2"),
    ("datatype", "i2"),
    ("bitpix", "i2"),
    ("slice_start", "i2"),
    ("pixdim", "f4", 8),
    ("vox_offset", "f4"),
    ("scl_slope", "f4"),
    ("scl_inter", "f4"),
    ("slice_end", "i2"),
    ("slice_code", "u1"),
    ("xyzt_units", "u1"),
    ("cal_max", "f4"),
    ("cal_min", "f4"),
    ("slice_duration", "f4"),
    ("toffset", "f4"),
    ("glmax", "i4"),
    ("glmin", "i4"),
    ("descrip", "S80"),
    ("aux_file", "S24"),
    ("qform_code", "i2"),
    ("sform_code", "i2"),
    ("quatern", "f4", 3),
    ("qoffset", "f4", 3),
    ("sform", "f4", [3, 4]),
    ("intent_name", "S16"),
    ("magic", "S4"),
])

HEADERTYPE2 = np.dtype([
    ("sizeof_hdr", "i4"),
    ("magic", "S8"),
    ("datatype", "i2"),
    ("bitpix", "i2"),
    ("dim", "i8", 8),
    ("intent_p", "f8", 3),
    ("pixdim", "f8", 8),
    ("vox_offset", "i8"),
    ("scl_slope", "f8"),
    ("scl_inter", "f8"),
    ("cal_max", "f8"),
    ("cal_min", "f8"),
    ("slice_duration", "f8"),
    ("toffset", "f8"),
    ("slice_start", "i8"),
    ("slice_end", "i8"),
    ("descrip", "S80"),
    ("aux_file", "S24"),
    ("qform_code", "i4"),
    ("sform_code", "i4"),
    ("quatern", "f8", 3),
    ("qoffset", "f8", 3),
    ("sform", "f8", [3, 4]),
    ("slice_code", "i4"),
    ("xyzt_units", "i4"),
    ("intent_code", "i4"),
    ("intent_name", "S16"),
    ("dim_info", "u1"),
    ("unused_str", "S15"),
])


class Recoder:

    def __init__(self, obj=None):
        self.forward = {}
        self.backward = {}
        if isinstance(obj, dict):
            for key, value in obj.items():
                self.forward[key] = value
                self.backward[value] = key
        elif isinstance(obj, (list, tuple)):
            for key, value in obj:
                self.forward[key] = value
                self.backward[value] = key

    def __getitem__(self, key):
        if key in self.forward:
            return self.forward[key]
        if key in self.backward:
            return self.backward[key]

    def __setitem__(self, key, value):
        self.forward[key] = value
        self.backward[value] = key

    def append(self, key_value: tuple):
        key, value = key_value
        self[key] = value

    def extend(self, key_values: list):
        for key, value in key_values:
            self[key] = value

    def update(self, key_values: dict):
        for key, value in key_values.items():
            self[key] = value


DTYPES = Recoder([
    (2, "uint8"),  # unsigned char (8 bits)
    (4, "int16"),  # signed short (16 bits)
    (8, "int32"),  # signed int (32 bits)
    (16, "single"),  # 32-bit float
    (32, "complex64"),  # 64-bit complex (2x 32-bit floats)
    (64, "double"),  # 64-bit float (double precision)
    (128, (("r", "uint8"), ("g", "uint8"), ("b", "uint8"))),  # 3x 8-bit unsigned char (RGB24)
    (256, "int8"),  # signed char (8 bits)
    (512, "uint16"),  # unsigned short (16 bits)
    (768, "uint32"),  # unsigned int (32 bits)
    (1024, "int64"),  # signed long long (64 bits)
    (1280, "uint64"),  # unsigned long long (64 bits)
    (1536, "double128"),  # 128-bit float (long double)
    (1792, "complex128"),  # 128-bit complex (2x 64-bit floats)
    (2048, "complex256"),  # 256-bit complex (2x 128-bit floats)
    (2304, (("r", "uint8"), ("g", "uint8"), ("b", "uint8"), ("a", "uint8"))),  # 4x 8-bit unsigned char (RGBA32)
])

UNITS = Recoder([
    (0, ""),
    (1, "m"),
    (2, "mm"),
    (3, "um"),
    (8, "s"),
    (16, "ms"),
    (24, "us"),
    (32, "hz"),
    (40, "ppm"),
    (48, "rad/s"),
])

JNIFTI_ZARR = Recoder([
    ("u1", "uint8"),
    ("i2", "int16"),
    ("i4", "int32"),
    ("f4", "single"),
    ("c8", "complex64"),
    ("f8", "double"),
    ((("r", "u1"), ("g", "u1"), ("b", "u1")), (("r", "uint8"), ("g", "uint8"), ("b", "uint8"))),
    ("i1", "int8"),
    ("u2", "uint16"),
    ("u4", "uint32"),
    ("i8", "int64"),
    ("u8", "uint64"),
    ("f16", "double128"),
    ("c16", "complex128"),
    ("c32", "complex256"),
    ((("r", "u1"), ("g", "u1"), ("b", "u1"), ("a", "u1")),
     (("r", "uint8"), ("g", "uint8"), ("b", "uint8"), ("a", "uint8"))),
    ("", ""),
    ("m", "meter"),
    ("mm", "millimeter"),
    ("um", "micrometer"),
    ("s", "second"),
    ("ms", "millisecond"),
    ("us", "microsecond"),
    ("hz", "hertz"),
    ("ppm", "micro"),
    ("rad/s", "radian"),
])

INTENTS = Recoder([
    # NIFTI_INTENT_NONE: Unknown data intent
    (0, ""),
    # NIFTI_INTENT_CORREL: Correlation coefficient R (1 param)
    (2, "corr"),
    # NIFTI_INTENT_TTEST: Student t statistic (1 param): p1 = DOF
    (3, "ttest"),
    # NIFTI_INTENT_FTEST: Fisher F statistic (2 params)
    (4, "ftest"),
    # NIFTI_INTENT_ZSCORE: Standard normal (0 params): Density = N(0,1)
    (5, "zscore"),
    # NIFTI_INTENT_CHISQ: Chi-squared (1 param): p1 = DOF
    (6, "chi2"),
    # NIFTI_INTENT_BETA: Beta distribution (2 params): p1=a, p2=b
    (7, "beta"),
    # NIFTI_INTENT_BINOM: Binomial distribution (2 params)
    (8, "binomial"),
    # NIFTI_INTENT_GAMMA: Gamma distribution (2 params)
    (9, "gamma"),
    # NIFTI_INTENT_POISSON: Poisson distribution (1 param): p1 = mean
    (10, "poisson"),
    # NIFTI_INTENT_NORMAL: Normal distribution (2 params)
    (11, "normal"),
    # NIFTI_INTENT_FTEST_NONC: Noncentral F statistic (3 params)
    (12, "ncftest"),
    # NIFTI_INTENT_CHISQ_NONC: Noncentral chi-squared statistic (2 params)
    (13, "ncchi2"),
    # NIFTI_INTENT_LOGISTIC: Logistic distribution (2 params)
    (14, "logistic"),
    # NIFTI_INTENT_LAPLACE: Laplace distribution (2 params)
    (15, "laplace"),
    # NIFTI_INTENT_UNIFORM: Uniform distribution: p1=lower end, p2=upper end
    (16, "uniform"),
    # NIFTI_INTENT_TTEST_NONC: Noncentral t statistic (2 params)
    (17, "ncttest"),
    # NIFTI_INTENT_WEIBULL: Weibull distribution (3 params)
    (18, "weibull"),
    # NIFTI_INTENT_CHI:Chi distribution (1 param): p1 = DOF
    (19, "chi"),
    # NIFTI_INTENT_INVGAUSS: Inverse Gaussian (2 params)
    (20, "invgauss"),
    # NIFTI_INTENT_EXTVAL: Extreme value type I (2 params)
    (21, "extval"),
    # NIFTI_INTENT_PVAL: Data is a 'p-value' (no params)
    (22, "pvalue"),
    # NIFTI_INTENT_LOGPVAL: Data is ln(p-value) (no params)
    (23, "logpvalue"),
    # NIFTI_INTENT_LOG10PVAL: Data is log10(p-value) (no params)
    (24, "log10pvalue"),
    # NIFTI_INTENT_ESTIMATE: Data is an estimate of some parameter
    (1001, "estimate"),
    # NIFTI_INTENT_LABEL: Data is an index into a set of labels
    (1002, "label"),
    # NIFTI_INTENT_NEURONAME: Data is an index into the NeuroNames labels
    (1003, "neuronames"),
    # NIFTI_INTENT_GENMATRIX: To store an M x N matrix at each voxel
    (1004, "matrix"),
    # NIFTI_INTENT_SYMMATRIX: To store an NxN symmetric matrix at each voxel
    (1005, "symmatrix"),
    # NIFTI_INTENT_DISPVECT: Each voxel is a displacement vector
    (1006, "dispvec"),
    # NIFTI_INTENT_VECTOR: Specifically for displacements
    (1007, "vector"),
    # NIFTI_INTENT_POINTSET: Any other type of vector
    (1008, "point"),
    # NIFTI_INTENT_TRIANGLE: Each voxel is really a triangle
    (1009, "triangle"),
    # NIFTI_INTENT_QUATERNION: Each voxel is a quaternion
    (1010, "quaternion"),
    # NIFTI_INTENT_DIMLESS: Dimensionless value - no params
    (1011, "unitless"),
    # NIFTI_INTENT_TIME_SERIES: Each data point is a time series
    (2001, "tseries"),
    # NIFTI_INTENT_NODE_INDEX: Each data point is a node index
    (2002, "elem"),
    # NIFTI_INTENT_RGB_VECTOR: Each data point is an RGB triplet
    (2003, "rgb"),
    # NIFTI_INTENT_RGBA_VECTOR: Each data point is a 4 valued RGBA
    (2004, "rgba"),
    # NIFTI_INTENT_SHAPE: Each data point is a shape value
    (2005, "shape"),
    (2006, "fsl_fnirt_displacement_field"),
    (2007, "fsl_cubic_spline_coefficients"),
    (2008, "fsl_dct_coefficients"),
    (2009, "fsl_quadratic_spline_coefficients"),
    (2016, "fsl_topup_cubic_spline_coefficients"),
    (2017, "fsl_topup_quadratic_spline_coefficients"),
    (2018, "fsl_topup_field"),
])

# parameters count for each intent
INTENTS_P = {key: 0 for key in INTENTS.backward.keys()}
INTENTS_P.update({
    "corr": 1,
    "ttest": 1,
    "ftest": 2,
    "chi2": 1,
    "beta": 2,
    "binomial": 2,
    "gamma": 2,
    "poisson": 1,
    "normal": 2,
    "ncftest": 3,
    "ncchi2": 2,
    "logistic": 2,
    "laplace": 2,
    "uniform": 2,
    "ncttest": 2,
    "weibull": 2,
    "chi": 1,
    "invgauss": 2,
    "extval": 2,
    "matrix": 2,
    "symmatrix": 1,
})

XFORMS = Recoder([
    (0, ""),
    (1, "scanner_anat"),
    (2, "aligned_anat"),
    (3, "talairach"),
    (4, "mni_152"),
    (5, "template_other"),
])

SLICEORDERS = Recoder([
    # NIFTI_SLICE_UNKNOWN: Unknown slice type, no parameters
    (0, ""),
    # NIFTI_SLICE_SEQ_INC: Slice sequentially increasing
    (1, "seq+"),
    # NIFTI_SLICE_SEQ_DEC: Slice sequentially decreasing
    (2, "seq-"),
    # NIFTI_SLICE_ALT_INC: Slice alternating increasing
    (3, "alt+"),
    # NIFTI_SLICE_ALT_DEC: Slice alternating decreasing
    (4, "alt-"),
    # NIFTI_SLICE_ALT_INC2: Slice alternating increasing type 2
    (5, "alt2+"),
    # NIFTI_SLICE_ALT_DEC2: Slice alternating decreasing type 2
    (6, "alt2-")
])


def get_magic_string(header):
    return re.sub(r'[\x00-\x1f]+', '', header['magic'].decode())


def try_header_version(buffer, version=1):
    if version == 1:
        HEADER_TYPE = HEADERTYPE1
        HEADER_SIZE = NIFTI_1_HEADER_SIZE
    elif version == 2:
        HEADER_TYPE = HEADERTYPE2
        HEADER_SIZE = NIFTI_2_HEADER_SIZE
    else:
        raise ValueError(f"Unsupported Nifti version {version}")
    byteorder_swapped = False
    header = np.frombuffer(buffer, dtype=HEADER_TYPE, count=1)[0]
    if header['sizeof_hdr'] != HEADER_SIZE:
        header = header.view(header.dtype.newbyteorder())
        byteorder_swapped = True
        if header['sizeof_hdr'] != HEADER_SIZE:
            return None

    validate_magic(header, version)
    return header, byteorder_swapped


def bin2nii(buffer, check_swapped=False):
    """
    Parameters
    ----------
    buffer : binary header data
    check_swapped : bool
        if true, return if the byte order is swapped
    Returns
    -------
    header : structure array of header
    swapped : bool

    """
    for v in (1, 2):
        result = try_header_version(buffer, v)
        if result:
            header, byteorder_swapped = result
            if check_swapped:
                return header, byteorder_swapped
            return header
    raise ValueError('Is this a nifti header?')


def validate_magic(header, version):
    magic_string = get_magic_string(header)
    if magic_string not in (f"n+{version}", f"ni{version}"):
        warnings.warn(f"Magic String {magic_string} does not match NIFTI version {version}")


# def get_header_version(buffer):
#     for v in (1, 2):
#         header = try_header_version(buffer, v)
#         if header:
#             return v
#     raise ValueError('Is this a nifti header?')
#
# def bin2nib(buffer):
#     version = get_header_version(buffer)
#     if version == 1:
#         NiftiHeader = Nifti1Header
#         NiftiImage = Nifti1Image
#     elif version==2:
#         NiftiHeader = Nifti2Header
#         NiftiImage = Nifti2Image
#     else:
#         raise ValueError(f"Unsupported Nifti version {version}")


def get_nibabel_klass(header):
    if header['sizeof_hdr'] == NIFTI_1_HEADER_SIZE:
        return Nifti1Header, Nifti1Image
    elif header['sizeof_hdr'] == NIFTI_2_HEADER_SIZE:
        return Nifti2Header, Nifti2Image
    else:
        raise ValueError(f"sizeof_hdr {header['sizeof_hdr']} does not match any Nifti header specification")
