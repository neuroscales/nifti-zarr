# NIfTI-Zarr draft specification

* Status of this document: release candidate
* Editor: Yael Balbastre <y.balbastre at ucl.ac.uk>
* Version: 1.0.rc1

## Abstract

This document specifies the nifti-zarr format for storing neuroimaging
data in the cloud.

## Table of Content

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
* [__JNIfTI__](https://github.com/NeuroJSON/jnifti/tree/master) is a pure JSON
  implementation of the NIfTI format.
* [__BIDS__](https://bids-specification.readthedocs.io) (Brain Imaging
  Data Structure) is a simple and intuitive way to organize and describe data.

## 1. Introduction

As biomedical imaging scales up, it is increasingly making use of remote
storage, remote computing and remote visualization. Classical file
formats—which store array data contiguously in a single file—are limited
at large scales as

1. They often do not store data at multiple resolutions;
2. They do not offer efficient parallel access to data chunks;
3. They do not efficiently compress 3D raster data.

These limits are very clear when it comes to visualizing very large data
volumes, which cannot be loaded in memory in full. In this context, it is
preferable to only load the data required to display a given scene (either
a large field-of-view at low-resolution, or a small field-of-view at high
resolution).

The [Zarr](https://zarr.readthedocs.io) format was developed to bypass
the limitations of single-file formats such as [HDF5](https://www.hdfgroup.org/).
The microscopy community is currently developping its own standard for
cloud-friendly biomedical imaging data
([OME-NGFF](https://ngff.openmicroscopy.org))—with a Zarr-based implementation—
and adds rules for storing multi-resolutions images and medical-specific
metadata such as axis names and voxel sizes. However, this community
has needs in terms of metadata and coordinate-space description that are
relatively complex, as they need to conform to different organs, a wide range
of acquisition systems, and different tissue processing pipelines. This has
drastically slowed down the adoption of a coordinate transform standard,
which hampers the use of OME-NGFF with neuroimaging data in two ways:

1. at the time of this writing, the current version of the format
   ([0.5](https://ngff.openmicroscopy.org/0.5/index.html)) only handles
   canonical scales and offsets;
2. the coordinate transform standard
   [being drafted](https://github.com/ome/ngff/pull/138) is more flexible
   than required for pure neuroimaging applications, which may prevent its
   widespread adoption by the neuroimaging community.

In contrast, the neuroimaging community has adopted and used a standard
"world" coordinate frame for decades, where

```text
+x = left      -> right
+y = posterior -> anterior
+z = inferior  -> superior
```

An affine transform is used to map from  the F-ordered voxel space $(i, j, k)$
to world space $(x, y, z)$. The neuroimaging community has also created a
simple data exchange format—[NIfTI](https://nifti.nimh.nih.gov)—that is
widely embraced and is the mandatory file format in standardization
efforts such as [BIDS](https://bids-specification.readthedocs.io).
However, the lack of multiresolution and/or chunk support in NIfTI has
led BIDS to adopt OME-TIFF and OME-ZARR as mandatory formats for its
microscopy component.

Let us further add that OME-NGFF is mostly oriented towards microscopy,
whereas neuroimaging formats focus on magnetic resonance imaging (MRI),
computed tomography (CT), and/or positon emission tomography (PET).

The NIfTI-Zarr (`nii.zarr`) specification attempts to merge the best of
both worlds, in the simplest possible way. Like NIfTI, it aims to make the
implementation of I/O libraries as simple as possible, to maximize chances
that it gets adopted by the community. Its guiding principles are

* __OME-Zarr compliant:__ any `nii.zarr` file should be a valid `ome.zarr` file;
* __OME-Zarr minimal:__ only implements the minimum set of metadata necessary
  to describe [multi-resolution] neuroimaging data;
* __NIfTI-compliant:__ the binary nifti header should be stored in its raw
  form in an additional `nifti` zarr array;
* __NIfTI-priority:__ if metadata conflict across the nifti header and OME
  attributes, the nifti metadata should take precedence.

> [!NOTE]
>
> * Being OME-NGFF compliant does not mean (for now) that the OME-NGFF
>   transform and the NIfTI transform match. Currently, OME-NGFF only handles
>   scales (for voxel sizes) and translations (for origin shifts caused by
>   pyramid methods). It is therefore impossible to encode an affine tranform -
>   or even swap axes - using the current OME-NGFF specification. What we
>   mean by OME-NGFF compliant is that any OME-NGFF viewer will correctly
>   display the content of the file _in scaled voxel space_.
> * OME-NGFF does not currently offer the possibility to store an intensity
>   transform. This means that OME-NGFF viewers will not use the intensity
>   affine transform encoded by `scl_slope` and `scl_inter` in the nifti header.
> * That said, the nifti layer added on top of OME-NGFF is light enough that
>   viewer developers may easily extend their software to handle
>   1. an affine geometric tranform, and
>   2. an affine intensity transform.
> * In modern languages such as Python and Julia, a virtual array
>   that points to the raw data can easily be encapsulated in a high-level
>   class that applies the intensity transform on the fly. This is
>   examplified in our Python and Julia reference implementations, which
>   respectively leverage `nibabel`'s `Nifti1Image` and `NIfTI.jl`'s `NIVolume`.

The simplicity of these guiding principles should make the adoption of
`nii.zarr` in cloud environments (almost) as straightforward as the adoption
of compressed-NIfTI (`.nii.gz`).

## 2. Format specification

A NIfTI-Zarr file __MUST__ be a valid
[OME-Zarr multi-resolution image](https://ngff.openmicroscopy.org/latest/#image-layout)
(and therefore also a valid
[Zarr dataset](https://zarr-specs.readthedocs.io/en/latest/specs.html)),
whose directory structure and metadata are described in sections
[2.1](#21-directory-structure), [2.2](#22-multiresolution-metadata) and
[2.3](#23-ome-ngff-metadata).

In addition, it __MUST__ store the nifti header corresponding to the finest
level of the pyramid as a Zarr array with the `"nifti"` key, as described
in section [2.4](#24-nifti-header).

> [!NOTE]
> At the time of this writing, Zarr is at version 3 and OME-NGFF is at
> version 0.5. NIfTI-Zarr being a simple layer within an OME-Zarr object,
> it does not impose a specific Zarr or OME-NGFF version. However, note that
> OME-NGFF v0.4 specifies that it should be used with Zarr v2, while
> OME-NGFF v0.5 specifies that it should be used with Zarr v3.
>
> Similarly, the NIfTI header may follow the
> [NIfTI v1](https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h) or
> [NIfTI v2](https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h)
> specifications.
>
> All examples below use the OME-NGFF v0.4 + Zarr v2 specification, but
> they could equivalently have used OME-NGFF v0.5 + Zarr v3.

### 2.1. Directory structure

__REF:__
- OME-NGFF latest: [https://ngff.openmicroscopy.org/latest/#image-layout](https://ngff.openmicroscopy.org/latest/#image-layout)
- OME-NGFF v0.5: [https://ngff.openmicroscopy.org/0.5/#image-layout](https://ngff.openmicroscopy.org/0.5/#image-layout)
- OME-NGFF v0.4: [https://ngff.openmicroscopy.org/0.4/#image-layout](https://ngff.openmicroscopy.org/0.4/#image-layout)

```text
└── mri.nii.zarr              # A NIfTI volume converted to Zarr.
    │
    ├── .zgroup               # Each volume is a Zarr group, of arrays.
    ├── .zattrs               # Group level attributes are stored in the .zattrs
    |                         # file include the OME "multiscales" key.
    |                         # In addition, the group level attributes
    │                         # may also contain "_ARRAY_DIMENSIONS" for
    |                         # compatibility with xarray if this group directly
    |                         # contains multi-scale arrays.
    |
    ├── nifti                 # The NIfTI header is stored as a raw array
    │   ├── .zarray           # of bytes, with an optional JSON form stored in
    │   ├── .zattrs           # the associated .zattrs file.
    │   └── 0
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

### 2.2. Multiresolution metadata

__REF:__
- Zarr v3: [https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html#array-metadata](https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html#array-metadata)
- Zarr v2: [https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html#arrays](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html#arrays)

```yaml
# {filename}.nii.zarr/{0..n}/.zarray
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
    "dtype": "<f4",         # SHOULD be the compatible with zattrs["nifti"]["DataType"]
    "fill_value": "NaN",    # Value to use for missing chunks
    "order": "F",           # MUST be "F"
    "shape": [              # SHOULD be the same as zattrs["nifti"]["Dim"][[3, 4, 2, 1, 0]]
        1,                  # T shape
        3,                  # C shape
        10000,              # Z shape
        10000,              # Y shape
        10000               # X shape
    ],
    "zarr_format": 2        # MUST be 2
}
```

### 2.3. OME-NGFF metadata

__REF:__
- OME-NGFF latest: [https://ngff.openmicroscopy.org/latest/#image-layout](https://ngff.openmicroscopy.org/latest/#multiscale-md)
- OME-NGFF v0.5: [https://ngff.openmicroscopy.org/0.5/#image-layout](https://ngff.openmicroscopy.org/0.5/#multiscale-md)
- OME-NGFF v0.4: [https://ngff.openmicroscopy.org/0.4/#image-layout](https://ngff.openmicroscopy.org/0.4/#multiscale-md)

```yaml
# {filename}.nii.zarr/.zattrs
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

__FUTURE CHANGES:__ as soon as affine coordinate transforms are integrated
into the OME-NGFF standard, it will be used to encode both the NIfTI
qform and sform as valid OME-NGFF metadata. However, only these two
transforms will be accepted in a valid NIfTI-Zarr file. None of the
more advanced combinations of affine and nonlinear transforms will be
accepted. Similarly, only a very specific set of coordinate spaces will
be accepted.

### 2.4. NIfTI header

__REF__:

* NIfTI v1: [https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h](https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h)
* NIfTI v2: [https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h](https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h)
* JNIfTI v1: [https://github.com/NeuroJSON/jnifti/blob/master/JNIfTI_specification.md#niftiheader](https://github.com/NeuroJSON/jnifti/blob/master/JNIfTI_specification.md#niftiheader)
* Zarr v3: [https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html#array-metadata](https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html#array-metadata)
* Zarr v2: [https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html#arrays](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html#arrays)
* NIfTI-Zarr v1.0.rc1: [nifti-zarr-schema-1.0.rc1.json](./nifti-zarr-schema-1.0.rc1.json)

The nifti header ([v1](https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h)
or [v2](https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h)) __MUST__ be encoded
as a single-chunk zarr array if bytes under the `nifti` key:

```text
└── mri.nii.zarr
    │
    └── nifti                 # The NIfTI header is stored as a raw array
        ├── .zarray           # of bytes, with an optional JSON form stored in
        ├── .zattrs           # the associated .zattrs file.
        └── 0
```

```yaml
# {filename}.nii.zarr/nifti/.zarray
{
    "zarr_format": 2,       # MUST be 2
    "shape": [348],         # MUST be header length if `dtype="u1"` else 1
    "compressor": {
        "id": "zlib",       # MUST be `{"id": "zlib", "level": 0..9}` or `null`
        "level": 9          # MUST be in 0..9, if "zlib"
    },
    "dtype": "u1",          # MUST be "u1" or "S{length:d}"
    "order": "C",           # UNUSED ("F" or "C")
    "chunks": [1],          # MUST be 1
    "fill_value": null,     # UNUSED
}
```

A JSON variant  of the nifti header __MAY__ be encoded in the attributes file
(`.zattrs` or `attrs.json`). The JSON variant is only provided for
human-readability. Its values __SHOULD__ be compatible with those of the
binary header. If values conflict between the binary and JSON headers,
the binary form __MUST__ take precedence. The JSON variant __MUST__ be
compatible with the [`JNIfTI/NIFTIHeader`](https://github.com/NeuroJSON/jnifti/blob/master/JNIfTI_specification.md#niftiheader)
specification. The JSON variant __MUST__ be compatible with the 
[`NIfTI-Zarr schema`](./nifti-zarr-schema-1.0.rc1.json) specification.

```yaml
# filename.nii.zarr/nifti/.zattrs
{
    # ------------------------------------------------------------------
    # All other tags **MAY** contain JSON representations of the nifti
    # header. Not all fields are included. This JSON representation is
    # OPTIONAL and has a lower priority than the binary header.
    # ------------------------------------------------------------------

    "NIIFormat": b"ni1\0",          # MUST be one of "ni1\0", "n+1\0", "ni2\0", "n+2\0"
    "Dim": [128, 128, 128, 1, 3],   # SHOULD match filename.nii.zarr/0/.zarray["shape"][[4, 3, 2, 0, 1]]
    "VoxelSize": [                  # xtztc unit size, SHOULD be compatible with:
        1.5, 1.5, 1.5,              #   filename.nii.zarr/.zattrs["multiscales"][0]["datasets"][0]["coordinateTransformations"][0]["scale"][2:5][::-1]
        0.1,                        #   filename.nii.zarr/.zattrs["multiscales"][0]["coordinateTransformations"][0]["scale"][0]
        1.0,                        #   filename.nii.zarr/.zattrs["multiscales"][0]["coordinateTransformations"][0]["scale"][1]
    ],
    "Unit": {                       # xyzt unit, SHOULD be compatible with:
        "L": "mm",                  #   filename.nii.zarr/.zattrs["multiscales"][0]["axes"][2:]["unit"]
        "T": "s",                   #   filename.nii.zarr/.zattrs["multiscales"][0]["axes"][0]["unit"]
    },
    "DataType": "single",           # SHOULD be compatible with filename.nii.zarr/0/.zarray["dtype"]
    "DimInfo": {
        "Freq": 1,                  # MUST be one of {0, 1, 2, 3}
        "Phase": 2,                 # MUST be one of {0, 1, 2, 3}
        "Slice": 3                  # MUST be one of {0, 1, 2, 3}
    },
    "Intent": "dispvec",            # MUST be a valid intent code (see table 4.2)
    "Param1": None,
    "Param2": None,
    "Param3": None,
    "Name": "",
    "ScaleSlope": 1.0,              # Data scaling: slope
    "ScaleOffset": 0.0,             # Data scaling: intercept
    "SliceType": "seq+",            # MUST be a valid slice timing code (see table 4.6)
    "FirstSliceID": 0 ,             # First slice index
    "LastSliceID": 127,             # Last slice index
    "SliceTime": 1.0,               # Time for 1 slice.
    "MinIntensity": 0.0,            # Min display intensity
    "MaxIntensity": 1.0,            # Max display intensity
    "TimeOffset": 0.0,              # Time axis shift
    "Description": "An MRI",        # Any text you like
    "AuxFile": "/path/to/aux",      # Auxiliary filename
    "QForm": 0,                     # MUST be a valid xform code (see table 4.5)
    "Quatern": {                    # Quaternion
        "b": b,
        "c": c,
        "d": d
    },
    "QuaternOffset": {              # Translation
        "x": tx,
        "y": ty,
        "z": tz
    },
    "SForm": 0,                     # MUST be a valid xform code (see table 4.5)
    "Affine": [
        [axx, axy, axz, tx],        # 1st row affine transform
        [ayx, ayy, ayz, ty],        # 2nd row affine transform
        [azx, azy, azz, tz],        # 3rd row affine transform
    ]
}
```

Some fields __SHOULD__ be equivalent to their OME-Zarr counterparts:

```python
nifti.zattrs["Dim"]            ==  0.zarray["shape"][[4, 3, 2, 0, 1]]   # Level 0 zarray
nifti.zattrs["DataType"]       ==  *.zarray["dtype"]                    # All zarrays
nifti.zattrs["VoxelSize"][:3]  ==  zattrs["multiscales"][0]["datasets"][0]["coordinateTransformations"][0]["scale"][2:5][::-1]
nifti.zattrs["VoxelSize"][3]   ==  zattrs["multiscales"][0]["coordinateTransformations"][0]["scale"][0]
nifti.zattrs["VoxelSize"][4]   ==  zattrs["multiscales"][0]["coordinateTransformations"][0]["scale"][1]
```

## 3. Main differences with NIfTI and/or OME-NGFF

* Following the OME-NGFF specifcation, dimensions are ordered as
  [T, C, Z, Y, X] (in C order) as opposed to [C, T, Z, Y, X].
* To conform with the NIfTI expectation, on load data should be returned
  as a [C, T, Z, Y, X] array (in C order; [X, Y, Z, T, C] in F order).

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

| Type       | Name             | NIfTI-1 usage                       | JNIfTI                | JSON type                              |
| ---------- | ---------------- | ----------------------------------- | --------------------- | -------------------------------------- |
| `int`      | `sizeof_hdr`     | __MUST__ be 348                     | `"NIIHeaderSize"`     | `integer`                              |
| `char`     | `data_type`      | ~~UNUSED~~                          | `"A75DataTypeName"`   | `string`                               |
| `char`     | `db_name`        | ~~UNUSED~~                          | `"A75DBName"`         | `string`                               |
| `int`      | `extents`        | ~~UNUSED~~                          | `"A75Extends"`        | `integer`                              |
| `short`    | `session_error`  | ~~UNUSED~~                          | `"A75SessionError"`   | `integer`                              |
| `char`     | `regular`        | ~~UNUSED~~                          | `"A75Regular"`        | `integer`                              |
| `char`     | `dim_info`       | MRI slice ordering.                 | `"DimInfo"`           | `property`                             |
|            |                  |                                     | `"DimInfo"/"Freq"`    | `integer`                              |
|            |                  |                                     | `"DimInfo"/"Phase"`   | `integer`                              |
|            |                  |                                     | `"DimInfo"/"Slice"`   | `integer`                              |
| `short[8]` | `dim`            | Data array dimensions.              | `"Dim"`               | `array[integer]`                       |
| `float`    | `intent_p1`      | 1st intent parameter.               | `"Param1"`            | `number`                               |
| `float`    | `intent_p2`      | 2nd intent parameter.               | `"Param2"`            | `number`                               |
| `float`    | `intent_p3`      | 3rd intent parameter.               | `"Param3"`            | `number`                               |
| `short`    | `intent_code`    | `NIFTI_INTENT_*` code.              | `"Intent"`            | `[integer, string]`                    |
| `short`    | `datatype`       | Defines data type!                  | `"DataType"`          | `[integer, string]`                    |
| `short`    | `bitpix`         | Number bits/voxel.                  | `"BitDepth"`          | `integer`                              |
| `short`    | `slice_start`    | First slice index.                  | `"FirstSliceID"`      | `integer`                              |
| `float[8]` | `pixdim`         | Grid spacings.                      | `"VoxelSize"`         | `array[number]`                        |
|            |                  |                                     | `"Orientation"`       | `property`                             |
|            |                  |                                     | `"Orientation"/"x"`   | `enum: ["l", "r", "a", "p", "i", "s"]` |
|            |                  |                                     | `"Orientation"/"y"`   | `enum: ["l", "r", "a", "p", "i", "s"]` |
|            |                  |                                     | `"Orientation"/"z"`   | `enum: ["l", "r", "a", "p", "i", "s"]` |
| `float`    | `vox_offset`     | Offset into .nii file               | `"NIIByteOffset"`     | `number`                               |
| `float`    | `scl_slope`      | Data scaling: slope.                | `"ScaleSlope"`        | `number`                               |
| `float`    | `scl_inter`      | Data scaling: offset.               | `"ScaleOffset"`       | `number`                               |
| `short`    | `slice_end`      | Last slice index.                   | `"LastSliceID"`       | `integer`                              |
| `char`     | `slice_code`     | Slice timing order.                 | `"SliceType"`         | `[integer, string]`                    |
| `char`     | `xyzt_units`     | Units of `pixdim[1..4]`             | `"Unit"`              | `property`                             |
|            |                  |                                     | `"Unit"/"L"`          | `[integer, string]`                    |
|            |                  |                                     | `"Unit"/"T"`          | `[integer, string]`                    |
| `float`    | `cal_max`        | Max display intensity               | `"MaxIntensity"`      | `number`                               |
| `float`    | `cal_min`        | Min display intensity               | `"MinIntensity"`      | `number`                               |
| `float`    | `slice_duration` | Time for 1 slice.                   | `"SliceTime"`         | `number`                               |
| `float`    | `toffset`        | Time axis shift.                    | `"TimeOffset"`        | `number`                               |
| `int`      | `glmax`          | ~~UNUSED~~                          | `"A75GlobalMax"`      | `integer`                              |
| `int`      | `glmin`          | ~~UNUSED~~                          | `"A75GlobalMin"`      | `integer`                              |
| `char[80]` | `descrip`        | any text you like.                  | `"Description"`       | `string`                               |
| `char[24]` | `aux_file`       | auxiliary filename.                 | `"AuxFile"`           | `string`                               |
| `short`    | `qform_code`     | `NIFTI_XFORM_*` code.               | `"QForm"`             | `integer`                              |
| `short`    | `sform_code`     | `NIFTI_XFORM_*` code.               | `"SForm"`             | `integer`                              |
| `float`    | `quatern_b`      | Quaternion b param.                 | `"Quatern"/"b"`       | `number`                               |
| `float`    | `quatern_c`      | Quaternion c param.                 | `"Quatern"/"c"`       | `number`                               |
| `float`    | `quatern_d`      | Quaternion d param.                 | `"Quatern"/"d"`       | `number`                               |
| `float`    | `qoffset_x`      | Quaternion x shift.                 | `"QuaternOffset"/"x"` | `number`                               |
| `float`    | `qoffset_y`      | Quaternion y shift.                 | `"QuaternOffset"/"y"` | `number`                               |
| `float`    | `qoffset_z`      | Quaternion z shift.                 | `"QuaternOffset"/"z"` | `number`                               |
| `float[4]` | `srow_x`         | 1st row affine transform.           | `"Affine"/[0]`        | `array[number]`                        |
| `float[4]` | `srow_y`         | 2nd row affine transform.           | `"Affine"/[1]`        | `array[number]`                        |
| `float[4]` | `srow_z`         | 3rd row affine transform.           | `"Affine"/[2]`        | `array[number]`                        |
| `char[16]` | `intent_name`    | 'name' or meaning of data.          | `"Name"`              | `string`                               |
| `char[4]`  | `magic`          | __MUST__ be `"ni1\0"` or `"n+1\0"`. | `"NIIFormat"`         | `string`                               |

### Table 4.2. Data types

In Zarr, the byte order __MUST__ be specified by prepending one of
`{"|", "<", ">"}` to the data type string.

| Data type    | NIfTI            | Zarr                                                                                      | JNIfTI         |
| ------------ | ---------------- | ----------------------------------------------------------------------------------------- | -------------- |
| `uint8`      | `2`              | <code>"&vert;u1"</code>                                                                   | `"uint8"`      |
| `int16`      | `4`              | `"i2"`                                                                                    | `"int16"`      |
| `int32`      | `8`              | `"i4"`                                                                                    | `"int32"`      |
| `float32`    | `16`             | `"f4"`                                                                                    | `"single"`     |
| `complex64`  | `32`             | `"c8"`                                                                                    | `"complex64"`  |
| `float64`    | `64`             | `"f8"`                                                                                    | `"double"`     |
| `rgb24`      | `128`            | <code>[["r", "&vert;u1"], ["g", "&vert;u1"], ["b", "&vert;u1"]]</code>                    | `"rgb24"`      |
| `int8`       | `256`            | <code>"&vert;i1"</code>                                                                   | `"int8"`       |
| `uint16`     | `512`            | `"u2"`                                                                                    | `"uint16"`     |
| `uint32`     | `768`            | `"u4"`                                                                                    | `"uint32"`     |
| `int64`      | `1024`           | `"i8"`                                                                                    | `"int64"`      |
| `uint64`     | `1280`           | `"u8"`                                                                                    | `"uint64"`     |
| `float128`   | `1536`           | `"f16"`                                                                                   | `"double128"`  |
| `complex128` | `1792`           | `"c16"`                                                                                   | `"complex128"` |
| `complex256` | `2048`           | `"c32"`                                                                                   | `"complex256"` |
| `rgba32`     | `2304`           | <code>[["r", "&vert;u1"], ["g", "&vert;u1"], ["b", "&vert;u1"], ["a", "&vert;u1"]]</code> | `"rgba32"`     |
| `bool`       | 🛑 unsupported! | <code>"&vert;b1"</code>                                                                   |                |
| `timedelta`  | 🛑 unsupported! | `"m8[{unit}]"`                                                                            |                |
| `time`       | 🛑 unsupported! | `"M8[{unit}]"`                                                                            |                |

### Table 4.3. Units

In OME-NGFF, units must be names from the UDUNITS-2 database.

| Unit        | NIfTI | UDUNITS-2       | OME-NGFF axis | JNIfTI    |
| ----------- | ----- | --------------- | ------------- | --------- |
| unknown     | `0`   | `""`            |               | `""`      |
| meter       | `1`   | `"meter"`       | `"space"`     | `"m"`     |
| millimeter  | `2`   | `"millimeter"`  | `"space"`     | `"mm"`    |
| micron      | `3`   | `"micrometer"`  | `"space"`     | `"um"`    |
| second      | `8`   | `"second"`      | `"time"`      | `"s"`     |
| millisecond | `16`  | `"millisecond"` | `"time"`      | `"ms"`    |
| microsecond | `24`  | `"microsecond"` | `"time"`      | `"us"`    |
| hertz       | `32`  | `"hertz"`       | `"channel"`   | `"hz"`    |
| ppm         | `40`  | `"micro"`       | `"channel"`   | `"ppm"`   |
| rad         | `48`  | `"radian"`      | `"channel"`   | `"rad/s"` |

### Table 4.4. Intents

| Intent                           | NIfTI  | JNIfTI          | `len(intent_p)` | Intent parameters                     |
| -------------------------------- | ------ | --------------- | --------------- | ------------------------------------- |
| None                             | `0`    | `"none"`        | `0`             | []                                    |
| Correlation coefficient R        | `2`    | `"corr"`        | `1`             | [dof]                                 |
| Student t statistic              | `3`    | `"ttest"`       | `1`             | [dof]                                 |
| Fisher F statistic               | `4`    | `"ftest"`       | `2`             | [num dof, den dof]                    |
| Standard normal                  | `5`    | `"zscore"`      | `0`             | []                                    |
| Chi-squared                      | `6`    | `"chi2"`        | `1`             | [dof]                                 |
| Beta distribution                | `7`    | `"beta"`        | `2`             | [a, b]                                |
| Binomial distribution            | `8`    | `"binomial"`    | `2`             | [nb trials, prob per trial]           |
| Gamma distribution               | `9`    | `"gamma"`       | `2`             | [shape, scale]                        |
| Poisson distribution             | `10`   | `"poisson"`     | `1`             | [mean]                                |
| Normal distribution              | `11`   | `"normal"`      | `2`             | [mean, standard deviation]            |
| Noncentral F statistic           | `12`   | `"ncftest"`     | `3`             | [num dof, den dof, num noncentrality] |
| Noncentral chi-squared statistic | `13`   | `"ncchi2"`      | `2`             | [dof, noncentrality]                  |
| Logistic distribution            | `14`   | `"logistic"`    | `2`             | [location, scale]                     |
| Laplace distribution             | `15`   | `"laplace"`     | `2`             | [location, scale]                     |
| Uniform distribution             | `16`   | `"uniform"`     | `2`             | [lower end, upper end]                |
| Noncentral t statistic           | `17`   | `"ncttest"`     | `2`             | [dof, noncentrality]                  |
| Weibull distribution             | `18`   | `"weibull"`     | `3`             | [location, scale, power]              |
| Chi distribution                 | `19`   | `"chi"`         | `1`             | [dof]                                 |
| Inverse Gaussian                 | `20`   | `"invgauss"`    | `2`             | [mu, lambda]                          |
| Extreme value type I             | `21`   | `"extval"`      | `2`             | [location, scale]                     |
| Data is a 'p-value'              | `22`   | `"pvalue"`      | `0`             | []                                    |
| Data is ln(p-value)              | `23`   | `"logpvalue"`   | `0`             | []                                    |
| Data is log10(p-value)           | `24`   | `"log10pvalue"` | `0`             | []                                    |
| Parameter estimate               | `1001` | `"estimate"`    | `0`             | []                                    |
| Index into set of labels         | `1002` | `"label"`       | `0`             | []                                    |
| Index into NeuroNames set        | `1003` | `"neuronames"`  | `0`             | []                                    |
| MxN matrix at each voxel         | `1004` | `"matrix"`      | `2`             | [M, N]                                |
| NxN matrix at each voxel         | `1005` | `"symmatrix"`   | `1`             | [N]                                   |
| Displacement field               | `1006` | `"dispvec"`     | `0`             | []                                    |
| Vector field                     | `1007` | `"vector"`      | `0`             | []                                    |
| Spatial coordinate               | `1008` | `"point"`       | `0`             | []                                    |
| Triangle (3 indices)             | `1009` | `"triangle"`    | `0`             | []                                    |
| Quaternion (4 values)            | `1010` | `"quaternion"`  | `0`             | []                                    |
| Dimensionless value              | `1011` | `"unitless"`    | `0`             | []                                    |
| Gifti time series                | `2001` | `"tseries"`     | `0`             | []                                    |
| Gifti node index                 | `2002` | `"elem"`        | `0`             | []                                    |
| Gifti RGB (3 values)             | `2003` | `"rgb"`         | `0`             | []                                    |
| Gifti RGBA (4 values)            | `2004` | `"rgba"`        | `0`             | []                                    |
| Gifti shape                      | `2005` | `"shape"`       | `0`             | []                                    |

Additional FSL codes:

| Intent                 | NIfTI  | JNIfTI                                      |
| ---------------------- | ------ | ------------------------------------------- |
| FSL displacement field | `2006` | `"FSL_FNIRT_DISPLACEMENT_FIELD"`            |
| FSL cubic spline       | `2007` | `"FSL_CUBIC_SPLINE_COEFFICIENTS"`           |
| FSL DCT coefficients   | `2008` | `"FSL_DCT_COEFFICIENTS"`                    |
| FSL quad spline        | `2009` | `"FSL_QUADRATIC_SPLINE_COEFFICIENTS"`       |
| FSL-TOPUP cubic spline | `2016` | `"FSL_TOPUP_CUBIC_SPLINE_COEFFICIENTS"`     |
| FSL-TOPUP quad spline  | `2017` | `"FSL_TOPUP_QUADRATIC_SPLINE_COEFFICIENTS"` |
| FSL-TOPUP field        | `2018` | `"FSL_TOPUP_FIELD"`                         |

### Table 4.5. Xforms

| Transform | NIfTI | JNIfTI        |
| --------- | ----- | ------------- |
| Unknown   | `0`   | `"unknown"`   |
| Scanner   | `1`   | `"scanner"`   |
| Aligned   | `2`   | `"aligned"`   |
| Talairach | `3`   | `"talairach"` |
| MNI       | `4`   | `"mni"`       |
| Template  | `5`   | `"template"`  |

### Table 4.6. Slice order

| Order                     | NIfTI | JNIfTI    |
| ------------------------- | ----- | --------- |
| Unknown                   | `0`   | `""`      |
| Sequential increasing     | `1`   | `"seq+"`  |
| Sequential decreasing     | `2`   | `"seq-"`  |
| alternating increasing    | `3`   | `"alt+"`  |
| alternating decreasing    | `4`   | `"alt-"`  |
| alternating increasing #2 | `5`   | `"alt2+"` |
| alternating decreasing #2 | `6`   | `"alt2-"` |

## 5. Reference implementations

Reference software to convert data between `.nii[.gz]` and `.nii.zarr`
is provided.

### 5.1. [Python](https://github.com/neuroscales/nifti-zarr-py)

```python
from niizarr import nii2zarr, zarr2nii

# convert from nii.gz to nii.zarr
nii2zarr('/path/to/mri.nii.gz', '/path/to/mri.nii.zarr')

# convert from nii.zarr to nii.gz
# The pyramid level can be selected with `level=L`, where 0 is the
# base/finest level.
zarr2nii('/path/to/mri.nii.zarr', '/path/to/mri.nii.gz')
zarr2nii('/path/to/mri.nii.zarr', '/path/to/mri_2x.nii.gz', level=1)

# Encapsulate a nifti-zarr into a nibabel.Nifti1Image object, whose
# dataobj is a dask array
img = zarr2nii('/path/to/mri.nii.zarr')
```

### 5.2. [Julia](https://github.com/neuroscales/NIfTIZarr.jl)

```julia
import NIfTIZarr: nii2zarr, zarr2nii

# convert from nii.gz to nii.zarr
nii2zarr("path/to/nifti.nii.gz", "s3://path/to/bucket")

# convert from nii.zarr to nii.gz
# The pyramid level can be selected with `level=L`, where 0 is the
# base/finest level.
zarr2nii("/path/to/mri.nii.zarr", "/path/to/mri.nii.gz")
zarr2nii("/path/to/mri.nii.zarr", "/path/to/mri_2x.nii.gz"; level=1)

# Encapsulate a nifti-zarr into a NIVolume object, whose
# raw data is a DiskArray
img = zarr2nii("s3://path/to/bucket")
```

### 5.3. [NGtools](https://github.com/neuroscales/ngtools)

The [`ngtools`](https://github.com/neuroscales/ngtools) package allows
NIfTI-Zarr files to be properly oriented and displayed in
[`neuroglancer`](https://github.com/google/neuroglancer).
