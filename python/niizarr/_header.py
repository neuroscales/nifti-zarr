import numpy as np


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

    def __init__(self, obj={}):
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
    (2, "u1"),
    (4, "i2"),
    (8, "i4"),
    (16, "f4"),
    (32, "c8"),
    (64, "f8"),
    (128, (("r", "u1"), ("g", "u1"), ("b", "u1"))),
    (256, "i1"),
    (512, "u2"),
    (768, "u4"),
    (1024, "i8"),
    (1280, "u8"),
    (1536, "f16"),
    (1792, "c16"),
    (2048, "c32"),
    (2304, (("r", "u1"), ("g", "u1"), ("b", "u1"), ("a", "u1"))),
])


UNITS = Recoder([
    (0, ""),
    (1, "meter"),
    (2, "millimeter"),
    (3, "micrometer"),
    (8, "second"),
    (16, "millisecond"),
    (24, "microsecond"),
    (32, "hertz"),
    (40, "micro"),
    (48, "radian"),
])

INTENTS = Recoder([
    (0, "NONE"),
    (2, "CORREL"),
    (3, "TTEST"),
    (4, "FTEST"),
    (5, "ZSCORE"),
    (6, "CHISQ"),
    (7, "BETA"),
    (8, "BINOM"),
    (9, "GAMMA"),
    (10, "POISSON"),
    (11, "NORMAL"),
    (12, "FTEST_NONC"),
    (13, "CHISQ_NONC"),
    (14, "LOGISTIC"),
    (15, "LAPLACE"),
    (16, "UNIFORM"),
    (17, "TTEST_NONC"),
    (18, "WEIBULL"),
    (19, "CHI"),
    (20, "INVGAUSS"),
    (21, "EXTVAL"),
    (22, "PVAL"),
    (23, "LOGPVAL"),
    (24, "LOG10PVAL"),
    (1001, "ESTIMATE"),
    (1002, "LABEL"),
    (1003, "NEURONAME"),
    (1004, "GENMATRIX"),
    (1005, "SYMMATRIX"),
    (1006, "DISPVECT"),
    (1007, "VECTOR"),
    (1008, "POINTSET"),
    (1009, "TRIANGLE"),
    (1010, "QUATERNION"),
    (1011, "DIMLESS"),
    (2001, "TIME_SERIES"),
    (2002, "NODE_INDEX"),
    (2003, "RGB_VECTOR"),
    (2004, "RGBA_VECTOR"),
    (2005, "SHAPE"),
    (2006, "FSL_FNIRT_DISPLACEMENT_FIELD"),
    (2007, "FSL_CUBIC_SPLINE_COEFFICIENTS"),
    (2008, "FSL_DCT_COEFFICIENTS"),
    (2009, "FSL_QUADRATIC_SPLINE_COEFFICIENTS"),
    (2016, "FSL_TOPUP_CUBIC_SPLINE_COEFFICIENTS"),
    (2017, "FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS"),
    (2018, "FSL_TOPUP_FIELD"),
])

INTENTS_P = {key: 0 for key in INTENTS.backward.keys()}
INTENTS_P.update({
    "CORREL": 1,
    "TTEST": 1,
    "FTEST": 2,
    "CHISQ": 1,
    "BETA": 2,
    "BINOM": 2,
    "GAMMA": 2,
    "POISSON": 1,
    "NORMAL": 2,
    "FTEST_NONC": 3,
    "CHISQ_NONC": 2,
    "LOGISTIC": 2,
    "LAPLACE": 2,
    "UNIFORM": 2,
    "TTEST_NONC": 2,
    "WEIBULL": 2,
    "CHI": 1,
    "INVGAUSS": 2,
    "EXTVAL": 2,
    "GENMATRIX": 2,
    "SYMMATRIX": 1,
})


XFORMS = Recoder([
    (0, "UNKNOWN"),
    (1, "SCANNER_ANAT"),
    (2, "ALIGNED_ANAT"),
    (3, "TALAIRACH"),
    (4, "MNI"),
    (5, "TEMPLATE_OTHER"),
])

SLICEORDERS = Recoder([
    (0, "UNKNOWN"),
    (1, "SEQ_INC"),
    (2, "SEQ_DEC"),
    (3, "ALT_INC"),
    (4, "ALT_DEC"),
    (5, "ALT_INC2"),
    (6, "ALT_DEC2"),
])
