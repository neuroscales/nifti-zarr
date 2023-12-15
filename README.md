# NIfTI-Zarr draft specification

## Abstract
This document contains _draft_ nifti-zarr specifications for storing
neuroimaging data in the cloud.

# Table of Content

0. [References](#0-references)
1. [Introduction](#1-introduction)
2. [Format Specification](#2-format-specification)
3. [Main differences with NIfTI and/or OME-NGFF](#3-main-differences-with-nifti-andor-ome-ngff)
4. [Conversion tables](#4-conversion-tables)
5. [Reference implementations](#5-reference-implementations)

## 0. References

* [__Zarr__](https://zarr.readthedocs.io) is a format for the storage of
  chunked, compressed, N-dimensional arrays inspired by HDF5, h5py and bcolz.
* [__OME-NGFF__](https://ngff.openmicroscopy.org) (Next Generation File
  Format) is a format based on zarr for the storage of biomedical imaging data.
* [__NIfTI__](https://nifti.nimh.nih.gov) (Neuroimaging Informatics
  Technology Initiative) is a single-file/single-resolution storage format
  for 3D+ neuroimaging data.
* [__BIDS__](https://bids-specification.readthedocs.io) (Brain Imaging
  Data Structure) is a simple and intuitive way to organize and describe data.


## 1. Introduction

As biomedical imaging scales up, it is making more and more use of remote
storage, remote computing and remote visualization. Classical file
formats—which store array data contiguously in a single file—are limited
at large scales as
1. They often do not store data at multiple resolutions;
2. They do not offer efficient parallel access to data chunks.
3. They do not efficiently compress 3D raster data

These limits are very clear when it comes to visualizing very large data
volumes, which cannot be loaded in full in memory. In this context, it is
preferable to only load the data required to display the current view (either
a large field-of-view at low-resolution, or a small field-of-view at high
resolution).

The zarr format was developed to bypass the limitations of single-file formats
such as HDF5. The microscopy community is currently developping its own
standard for cloud-friendly biomedical imaging data: OME-NGFF. It build on
zarr and adds rules for storing multi-resolutions images and medical-specific
metadata such as axis names and voxel sizes. However, the microscopy community
has needs in terms of metadata and coordinate-space description that are
relatively complex, as they need to conform to different organs, different
acquisition systems, or different tissue processing pipelines. This has
drastically slowed down the adoption of a coordinate transform standard—the
current version (0.5) only handles canonical scales and offsets—and has also
made the future coordinate transform standard much more complicated.

In contrast, the neuroimaging community has adopted and used a standard
"world" coordinate frame for decades, where
```
+x = left      -> right
+y = posterior -> anterior
+z = inferior  -> superior
```
An affine transform is used to map from  the F-ordered voxel space (i, j, k)
to world space (x, y, z). The neuroimaging community has also created a
simple data exchange format—NIfTI—that has been widely embraced and is the
mandatory file format in standardization efforts such as BIDS. However, the
lack of multiresolution/chunk support in NIfTI has lead BIDS to adopt OME-TIFF
and OME-ZARR as standard format for its microscopy component.

The NIfTI-Zarr (`nii.zarr`) specification attempts to merge the best of
both worlds, in the simplest possible way. Like NIfTI, it aims to make the
implementation of I/O libraries as simple as possible, to maximize chances
that it gets adopted by the community. Its guiding principles are
* __OME-NGFF compliant:__ any `nii.zarr` file should be a valid `ome.zarr` file.
* __OME-NGFF minimal:__ only implements the minimum set of metadata necessary
  to describe [multi-resolution] neuroimaging data
* __NIfTI-complicant:__ the binary nifti header should be stored in the
  [group-level] `.zattrs` object.
* __NIfTI-priority:__ if metadata conflict across the nifti header and OME
  attributes, the nifti metadata should take precedence.

**NOTES**
* Being OME-NGFF compliant does not mean (for now) that the OME-NGFF
  transform and the NIfTI transform match. Currently, OME-NGFF only handles
  scales (for voxel sizes) and translation (for origin shifts caused by
  pyramid methods). It is therefore impossible to encode an affine tranform -
  or even swap axes - using the current OME-NGFF specification. What we
  mean by OME-NGFF compliant is that any OME-NGFF viewer will correctly
  display the content of the file _in scaled voxel space_.
* OME-NGFF does not currently offer the possibility to store an intensity
  transform. This means that OME-NGFF viewers will not use the intensity
  affine transform encoded by `scl_slope` and `scl_inter` in the nifti header.
* That said, the nifti layer added on top of OME-NGFF is light enough that
  viewer developers may easily extend their software to handle
  1. an affine geometric tranform, and
  2. an affine intensity transform.

The simplicity of these guiding principles should make the adoption of
`nii.zarr` in cloud environments (almost) as straightforward as the adoption
of compressed-NIfTI (`.nii.gz`).

## 2. Format specification

### 2.1. Directory structure
Directory structure:
```
└── mri.nii.zarr              # A nifti volume converted to Zarr.
    │
    ├── .zgroup               # Each volume is a Zarr group, of arrays.
    ├── .zattrs               # Group level attributes are stored in the .zattrs
    |                         # file and include "multiscales" and "nifti" (see
    |                         # below). In addition, the group level attributes
    │                         # may also contain "_ARRAY_DIMENSIONS" for
    |                         # compatibility with xarray if this group directly
    |                         # contains multi-scale arrays.
    │
    ├── 0                     # Each multiscale level is stored as a separate
    |                         # Zarr array, which is a folder containing chunk
    │   ...                   # files which compose the array.
    |
    └── n                     # The name of the array is arbitrary with the
        |                     # ordering defined by the "multiscales" metadata,
        │                     # but is often a sequence starting at 0.
        │
        ├── .zarray           # All image arrays must be up to 5-dimensional
        │                     # with the axis of type time before type channel,
        |                     # before spatial axes.
        │
        └─ t                  # Chunks are stored with the nested directory
           └─ c               # layout. All but the last chunk element are stored
              └─ z            # as directories. The terminal chunk is a file.
                 └─ y         # Together the directory and file names provide the
                    └─ x      # "chunk coordinate" (t, c, z, y, x), where the
                              # maximum coordinate will be dimension_size / chunk_size.
```

### 2.2. Zarr metadata

**REF:** https://zarr.readthedocs.io/en/stable/spec/v2.html#arrays

`.zarray`
```python
{
    "chunks": [
        1,                  # Number of time chunks
        3,                  # Number of channel chunks
        1000,               # Number of z chunks
        1000,               # Number of y chunks
        1000,               # Number of x chunks
    ],
    "compressor": {         # Codec used to compress chunks
        "id": "blosc",      # MUST be "blosc" or "zlib"
        "cname": "lz4",
        "clevel": 5,
        "shuffle": 1
    },
    "dtype": "<f4",         # SHOULD be the same as zattrs["nifti"]["datatype"]
    "fill_value": "NaN",    # Value to use for missing chunks
    "order": "F",           # MUST be "F"
    "shape": [              # SHOULD be the same as zattrs["nifti"]["dim"][[3, 4, 2, 1, 0]]
        1,                  # T shape
        3,                  # C shape
        10000,              # Z shape
        10000,              # Y shape
        10000               # X shape
    ],
    "zarr_format": 2
}
```

### 2.3. OME-NGFF metadata

**REF:** https://ngff.openmicroscopy.org/0.4/#multiscale-md

`.zattrs`
```python
{
    "multiscales": [
        {
            "version": "0.4",
            "axes": [
                {"name": "t", "type": "time", "unit": "second"},
                {"name": "c", "type": "channel"},
                {"name": "z", "type": "space", "unit": "millimeter"},
                {"name": "y", "type": "space", "unit": "millimeter"},
                {"name": "x", "type": "space", "unit": "millimeter"}
            ],
            "datasets": [
                {
                    "path": "0",
                    "coordinateTransformations": [{
                        # the voxel size for the first scale level (0.5 millimeter)
                        "type": "scale",
                        "scale": [1.0, 1.0, 0.5, 0.5, 0.5]
                    }]
                },
                {
                    "path": "1",
                    "coordinateTransformations": [{
                        # the voxel size for the second scale level
                        # (downscaled by a factor of 2 -> 1 millimeter)
                        "type": "scale",
                        "scale": [1.0, 1.0, 1.0, 1.0, 1.0]
                    }]
                },
                {
                    "path": "2",
                    "coordinateTransformations": [{
                        # the voxel size for the third scale level
                        # (downscaled by a factor of 4 -> 2 millimeter)
                        "type": "scale",
                        "scale": [1.0, 1.0, 2.0, 2.0, 2.0]
                    }]
                }
            ],
            "coordinateTransformations": [{
                # the time unit (0.1 seconds),
                # which is the same for each scale level
                "type": "scale",
                "scale": [0.1, 1.0, 1.0, 1.0, 1.0]
            }],
        }
    ]
}
```

### 2.4. NIfTI header

**REF**: https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h

The nifti header ([v1](https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h)
or [v2](https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h)) **MUST** be encoded as
a string using [base64](https://en.wikipedia.org/wiki/Base64) encoding and saved in
the group-level `.zattrs` under the `["nifti"]["base64"]` key:

A JSON version  of the nifti header **MAY** be encoded under the `"nifti"` key.
The JSON version is only provided for human-readability. Its values **SHOULD** be
compatible with those of the binary header. If values conflict between the binary
and JSON headers, the binary form **MUST** take precedence.

`.zattrs`
```python
{
   "nifti": {
      "base64": "...",               # MUST be present. Base64 encoding of the binary header.

      # All other tags **MAY** contain JSON representations
      # of the nifti header. Not all fields are included.
      # This JSON representation is **optional** and has
      # a lower priority than the base64 header.

      "magic": b"nz1\0",             # MUST be "nz1\0" or "nz2\0"
      "dim" : [128, 128, 128, 1, 3], # SHOULD match .zarray["shape"][[4, 3, 2, 0, 1]]
      "pixdim": [                    # xtztc unit size, SHOULD match:
         1.5, 1.5, 1.5,              #   .zattrs["multiscales"][0]["datasets"][0]["coordinateTransformations"][-1]["scale"][2:5]
         0.1,                        #   .zattrs["multiscales"][0]["coordinateTransformations"][-1]["scale"][0]
         1.0,                        #   .zattrs["multiscales"][0]["coordinateTransformations"][-1]["scale"][1]
      ],
      "units": {                     # xyzt unit, SHOULD match
         "space": "millimeter",      #   .zattrs["multiscales"][0]["axes"][2:]["unit"]
         "time": "second",           #   .zattrs["multiscales"][0]["axes"][0]["unit"]
      },
      "datatype": "<f4",             # MUST be a zarr-compatible data type
                                     # SHOULD match .zarray["dtype"]
      "dim_info": {
         "freq": 1,                  # MUST be one of {0, 1, 2, 3}
         "phase": 2,                 # MUST be one of {0, 1, 2, 3}
         "slice": 3                  # MUST be one of {0, 1, 2, 3}
      },
      "intent": {
         "code": "DISPVECT",         # MUST be a valid intent code (see table 4.2)
         "name": "",                 # 'name' or meaning of data
         "p": []                     # Intent parameters (see table 4.2)
      },
      "scl": {
         "slope": 1.0,               # Data scaling: slope
         "inter": 0.0                # Data scaling: intercept
      },
      "slice": {                     # Slice timing order
         "code": "SEQ_INC",          # MUST be a valid slice timing code (see table 4.6)
         "start": 0 ,                # First slice index
         "end": 127,                 # Last slice index
         "duration": 1.0             # Time for 1 slice.
      },
      "cal": {
         "min": 0.0,                 # Min display intensity
         "max": 1.0                  # Max display intensity
      },
      "toffset": 0.0,                # Time axis shift
      "description": "An MRI",       # Any text you like
      "aux_file": "/path/to/aux",    # Auxiliary filename
      "qform": {
         "code": "SCANNER_ANAT",     # MUST be a valid xform name (see table 4.5)
         "quatern": [b, c, d],       # Quaternion
         "offset": [tx, ty, tz]      # Translation
      },
      "sform": {
         "code": "ALIGN_ANAT",       # MUST be a valid xform name (see table 4.5)
         "affine": [
            [axx, axy, axz, tx],     # 1st row affine transform
            [ayx, ayy, ayz, ty],     # 2nd row affine transform
            [azx, azy, azz, tz],     # 3rd row affine transform
         ]
      },
  }
}
```

Some fields **SHOULD** be equivalent to their OME-Zarr counterparts:
```python
zattrs["nifti"]["dim"]         ==  zarray["shape"][[4, 3, 2, 0, 1]]   # Level 0 zarray
zattrs["nifti"]["datatype"]    ==  zarray["dtype"]                    # All zarrays
zattrs["nifti"]["pixdim"][:3]  ==  zattrs["multiscales"][0]["datasets"][0]["coordinateTransformations"][-1]["scale"][2:5::-1]
zattrs["nifti"]["pixdim"][3]   ==  zattrs["multiscales"][0]["coordinateTransformations"][-1]["scale"][0]
zattrs["nifti"]["pixdim"][4]   ==  zattrs["multiscales"][0]["coordinateTransformations"][-1]["scale"][1]
```

## 3. Main differences with NIfTI and/or OME-NGFF

* Following the OME-NGFF specifcation, dimensions are ordered as [T, C, Z, Y, X] (in C order)
  as opposed to [C, T, Z, Y, X].
* To conform with the NIfTI expectation, on load data should be returned as a [C, T, Z, Y, X] array
  (in C order; [X, Y, Z, T, C] in F order).

### 3.1. NIfTI features that are not supported by NIfTI-Zarr

* Any file with more than 5 dimensions

### 3.2. OME-NGFF features that are not supported by NIfTI-Zarr

* Image collections
* Image with labels
* High-content screening data
* OMERO metadata

## 4. Conversion tables

### Table 4.1. NIfTI header

As a reminder, the nifti1 header has the following structure:
| Type       | Name             | NIFfI-1 usage                   |
| ---------- | ---------------- | ------------------------------- |
| `int`      | `sizeof_hdr`     | **MUST** be 348                 |
| `char`     | `data_type`      | ~~UNUSED~~                      |
| `char`     | `db_name`        | ~~UNUSED~~                      |
| `int`      | `extents`        | ~~UNUSED~~                      |
| `short`    | `session_error`  | ~~UNUSED~~                      |
| `char`     | `regular`        | ~~UNUSED~~                      |
| `char`     | `dim_info`       | MRI slice ordering.             |
| `short[8]` | `dim`            | Data array dimensions.          |
| `float`    | `intent_p1`      | 1st intent parameter.           |
| `float`    | `intent_p2`      | 2nd intent parameter.           |
| `float`    | `intent_p3`      | 3rd intent parameter.           |
| `short`    | `intent_code`    | `NIFTI_INTENT_*` code.          |
| `short`    | `datatype`       | Defines data type!              |
| `short`    | `bitpix`         | Number bits/voxel.              |
| `short`    | `slice_start`    | First slice index.              |
| `float[8]` | `pixdim`         | Grid spacings.                  |
| `float`    | `vox_offset`     | Offset into .nii file           |
| `float`    | `scl_slope`      | Data scaling: slope.            |
| `float`    | `scl_inter`      | Data scaling: offset.           |
| `short`    | `slice_end`      | Last slice index.               |
| `char`     | `slice_code`     | Slice timing order.             |
| `char`     | `xyzt_units`     | Units of `pixdim[1..4]`         |
| `float`    | `cal_max`        | Max display intensity           |
| `float`    | `cal_min`        | Min display intensity           |
| `float`    | `slice_duration` | Time for 1 slice.               |
| `float`    | `toffset`        | Time axis shift.                |
| `int`      | `glmax`          | ~~UNUSED~~                      |
| `int`      | `glmin`          | ~~UNUSED~~                      |
| `char[80]` | `descrip`        | any text you like.              |
| `char[24]` | `aux_file`       | auxiliary filename.             |
| `short`    | `qform_code`     | `NIFTI_XFORM_*` code.           |
| `short`    | `sform_code`     | `NIFTI_XFORM_*` code.           |
| `float`    | `quatern_b`      | Quaternion b param.             |
| `float`    | `quatern_c`      | Quaternion c param.             |
| `float`    | `quatern_d`      | Quaternion d param.             |
| `float`    | `qoffset_x`      | Quaternion x shift.             |
| `float`    | `qoffset_y`      | Quaternion y shift.             |
| `float`    | `qoffset_z`      | Quaternion z shift.             |
| `float[4]` | `srow_x`         | 1st row affine transform.       |
| `float[4]` | `srow_y`         | 2nd row affine transform.       |
| `float[4]` | `srow_z`         | 3rd row affine transform.       |
| `char[16]` | `intent_name`    | 'name' or meaning of data.      |
| `char[4]`  | `magic`          | **MUST** be `"ni1\0"` or `"n+1\0"`. |


### Table 4.2. Data types

In Zarr, the byte order **MUST** be specified by prepending one of `{"|", "<", ">"}` to the data type string.
| Data type    | NIfTI  | Zarr    |
| ------------ | ------ | ------- |
| `uint8`      | `2`    | <code>"&vert;u1"</code> |
| `int16`      | `4`    | `"i2"`  |
| `int32`      | `8`    | `"i4"`  |
| `float32`    | `16`   | `"f4"`  |
| `complex64`  | `32`   | `"c8"`  |
| `float64`    | `64`   | `"f8"`  |
| `rgb24`      | `128`  | <code>[["r", "&vert;u1"], ["g", "&vert;u1"], ["b", "&vert;u1"]]</code> |
| `int8`       | `256`  | <code>"&vert;i1"</code> |
| `uint16`     | `512`  | `"u2"`  |
| `uint32`     | `768`  | `"u4"`  |
| `int64`      | `1024` | `"i8"`  |
| `uint64`     | `1280` | `"u8"`  |
| `float128`   | `1536` | `"f16"` |
| `complex128` | `1792` | `"c16"` |
| `complex256` | `2048` | `"c32"` |
| `rgba32`     | `2304` | <code>[["r", "&vert;u1"], ["g", "&vert;u1"], ["b", "&vert;u1"], ["a", "&vert;u1"]]</code> |
| `bool`       | 🛑 unsupported! | <code>"&vert;b1"</code> |
| `timedelta`  | 🛑 unsupported! | `"m8[{unit}]"` |
| `time`       | 🛑 unsupported! | `"M8[{unit}]"` |

### Table 4.3. Units

In OME-NGFF, units must be names from the UDUNITS-2 database.
| Unit         | NIfTI  | UDUNITS-2       | OME-NGFF axis |
| ------------ | ------ | --------------- | ------------- |
| unknown      | `0`    | `""`            |               |
| meter        | `1`    | `"meter"`       | `"space"`     |
| millimeter   | `2`    | `"millimeter"`  | `"space"`     |
| micron       | `3`    | `"micrometer"`  | `"space"`     |
| second       | `8`    | `"second"`      | `"time"`      |
| millisecond  | `16`   | `"millisecond"` | `"time"`      |
| microsecond  | `24`   | `"microsecond"` | `"time"`      |
| hertz        | `32`   | `"hertz"`       | `"channel"`   |
| ppm          | `40`   | `"micro"`       | `"channel"`   |
| rad          | `48`   | `"radian"`      | `"channel"`   |

### Table 4.4. Intents

| Intent                    | NIfTI  | NIfTI-Zarr's JSON header | `len(intent_p)` | Intent parameters |
| ------------------------- | ------ | ------------------------ | --------------- | ----------------- |
| None                      | `0`    | `"NONE"`                 | `0` | [] |
| Correlation coefficient R | `2`    | `"CORREL"`               | `1` | [dof] |
| Student t statistic       | `3`    | `"TTEST"`                | `1` | [dof] |
| Fisher F statistic        | `4`    | `"FTEST"`                | `2` | [num dof, den dof] |
| Standard normal           | `5`    | `"ZSCORE"`               | `0` | [] |
| Chi-squared               | `6`    | `"CHISQ"`                | `1` | [dof] |
| Beta distribution         | `7`    | `"BETA"`                 | `2` | [a, b] |
| Binomial distribution     | `8`    | `"BINOM"`                | `2` | [nb trials, prob per trial] |
| Gamma distribution        | `9`    | `"GAMMA"`                | `2` | [shape, scale] |
| Poisson distribution      | `10`   | `"POISSON"`              | `1` | [mean] |
| Normal distribution       | `11`   | `"NORMAL"`               | `2` | [mean, standard deviation] |
| Noncentral F statistic    | `12`   | `"FTEST_NONC"`           | `3` | [num dof, den dof, num noncentrality] |
| Noncentral chi-squared statistic | `13`  | `"CHISQ_NONC"`     | `2` | [dof, noncentrality] |
| Logistic distribution     | `14`   | `"LOGISTIC"`             | `2` | [location, scale] |
| Laplace distribution      | `15`   | `"LAPLACE"`              | `2` | [location, scale] |
| Uniform distribution      | `16`   | `"UNIFORM"`              | `2` | [lower end, upper end] |
| Noncentral t statistic    | `17`   | `"TTEST_NONC"`           | `2` | [dof, noncentrality] |
| Weibull distribution      | `18`   | `"WEIBULL"`              | `3` | [location, scale, power] |
| Chi distribution          | `19`   | `"CHI"`                  | `1` | [dof] |
| Inverse Gaussian          | `20`   | `"INVGAUSS"`             | `2` | [mu, lambda] |
| Extreme value type I      | `21`   | `"EXTVAL"`               | `2` | [location, scale] |
| Data is a 'p-value'       | `22`   | `"PVAL"`                 | `0` | [] |
| Data is ln(p-value)       | `23`   | `"LOGPVAL"`              | `0` | [] |
| Data is log10(p-value)    | `24`   | `"LOG10PVAL"`            | `0` | [] |
| Parameter estimate        | `1001` | `"ESTIMATE"`             | `0` | [] |
| Index into set of labels  | `1002` | `"LABEL"`                | `0` | [] |
| Index into NeuroNames set | `1003` | `"NEURONAME"`            | `0` | [] |
| MxN matrix at each voxel  | `1004` | `"GENMATRIX"`            | `2` | [M, N] |
| NxN matrix at each voxel  | `1005` | `"SYMMATRIX"`            | `1` | [N] |
| Displacement field        | `1006` | `"DISPVECT"`             | `0` | [] |
| Vector field              | `1007` | `"VECTOR"`               | `0` | [] |
| Spatial coordinate        | `1008` | `"POINTSET"`             | `0` | [] |
| Triangle (3 indices)      | `1009` | `"TRIANGLE"`             | `0` | [] |
| Quaternion (4 values)     | `1010` | `"QUATERNION"`           | `0` | [] |
| Dimensionless value       | `1011` | `"DIMLESS"`              | `0` | [] |
| Gifti time series         | `2001` | `"TIME_SERIES"`          | `0` | [] |
| Gifti node index          | `2002` | `"NODE_INDEX"`           | `0` | [] |
| Gifti RGB (3 values)      | `2003` | `"RGB_VECTOR"`           | `0` | [] |
| Gifti RGBA (4 values)     | `2004` | `"RGBA_VECTOR"`          | `0` | [] |
| Gifti shape               | `2005` | `"SHAPE"`                | `0` | [] |

Additional FSL codes:
| Intent                    | NIfTI  | NIfTI-Zarr's JSON header                    |
| ------------------------- | ------ | ------------------------------------------- |
| FSL displacement field    | `2006` | `"FSL_FNIRT_DISPLACEMENT_FIELD"`            |
| FSL cubic spline          | `2007` | `"FSL_CUBIC_SPLINE_COEFFICIENTS"`           |
| FSL DCT coefficients      | `2008` | `"FSL_DCT_COEFFICIENTS"`                    |
| FSL quad spline           | `2009` | `"FSL_QUADRATIC_SPLINE_COEFFICIENTS"`       |
| FSL-TOPUP cubic spline    | `2016` | `"FSL_TOPUP_CUBIC_SPLINE_COEFFICIENTS"`     |
| FSL-TOPUP quad spline     | `2017` | `"FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS"` |
| FSL-TOPUP field           | `2018` | `"FSL_TOPUP_FIELD"`                         |

### Table 4.5. Xforms

| Transform | NIfTI | NIfTI-Zarr's JSON header |
| --------- | ----- | ------------------------ |
| Unknown   | `0`   | `"UNKNOWN"`              |
| Scanner   | `1`   | `"SCANNER_ANAT"`         |
| Aligned   | `2`   | `"ALIGNED_ANAT"`         |
| Talairach | `3`   | `"TALAIRACH"`            |
| MNI       | `4`   | `"MNI"`                  |
| Template  | `5`   | `"TEMPLATE_OTHER"`       |

### Table 4.6. Slice order

| Order                     | NIfTI | NIfTI-Zarr's JSON header |
| ------------------------- | ----- | ------------------------ |
| Unknown                   | `0`   | `"UNKNOWN"`              |
| Sequential increasing     | `1`   | `"SEQ_INC"`              |
| Sequential decreasing     | `2`   | `"SEQ_DEC"`              |
| alternating increasing    | `3`   | `"ALT_INC"`              |
| alternating decreasing    | `4`   | `"ALT_DEC"`              |
| alternating increasing #2 | `5`   | `"ALT_INC2"`             |
| alternating decreasing #2 | `6`   | `"ALT_DEC2"`             |


## 5. Reference implementations

We implemented software to convert data between `.nii[.gz]` and `.nii.zarr`.

### 5.1. Python

`python/niizarr`
```python
from niizarr import nii2zarr, zarr2nii

# convert from nii.gz to nii.zarr
nii2zarr('/path/to/mri.nii.gz', '/path/to/mri.nii.zarr')

# convert from nii.zarr to nii.gz
zarr2nii('/path/to/mri.nii.zarr', '/path/to/mri.nii.gz')
zarr2nii('/path/to/mri.nii.zarr', '/path/to/mri_2x.nii.gz', level=1)

# Return a nibabel.Nifti1Image or nibabel.Nifti2Image, whose dataobj is
# a dask array
img = zarr2nii('/path/to/mri.nii.zarr')
```

Example script loading a nii.zarr in a neuroglancer instance and
setting the correct affine.
```python
import neuroglancer as ng
import zarr
import nibabel
import niizarr
import io
import base64
import numpy as np
import fsspec
import json

URL = 'https://path/to/mri.nii.zarr'

with fsspec.open(URL + '/.zattrs', 'r') as f:
    zattrs = json.load(f)

# parse header with nibabel (we must fix magic string first)
hdr = niizarr.bin2nii(base64.b64decode(zattrs["nifti"]["base64"])).copy()
hdr["magic"] = "n+1"
hdr = nibabel.Nifti1Header.from_fileobj(io.BytesIO(hdr))

# compute matrices
permute = np.eye(4)
permute[:-1, :-1] = permute[:-1, [2, 1, 0]]
vox2world = hdr.get_best_affine()
vox2phys = np.eye(4)
scale = zattrs["multiscales"][0]["datasets"][0]["coordinateTransformations"][0]["scale"]
shift = zattrs["multiscales"][0]["datasets"][0]["coordinateTransformations"][-1]["translation"]
vox2phys[[0, 1, 2], [0, 1, 2]] = list(reversed(scale[2:]))
vox2phys[:-1, -1] = list(reversed(shift[2:]))
phys2world = vox2world @ permute @ np.linalg.inv(vox2phys)

# define neuroglancer transform
ras_space = ng.CoordinateSpace(
    names=["x", "y", "z"],
    units="mm",
    scales=[1]*3,
)
phys_space = ng.CoordinateSpace(
    names=["z", "y", "x"],
    units="mm",
    scales=scale[2:],
)
transform = ng.CoordinateSpaceTransform(
    matrix=phys2world[:3, :4],
    input_dimensions=phys_space,
    output_dimensions=ras_space,
)

# launch neuroglancer instance
viewer = ng.Viewer()
print(viewer.get_viewer_url())

# load volume and transform
with viewer.txn() as state:
    state.layers.append(
        name="mri",
        layer=ng.ImageLayer(
                source=ng.LayerDataSource(
                url="zarr://" + URL,
                transform=transform
            )
        )
    )
```
