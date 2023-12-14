# NIfTI-Zarr draft specification

## Abstract
This document contains _draft_ nifti-zarr specifications for storing neuroimaging data in the cloud.

## References

* [__Zarr__](https://zarr.readthedocs.io) is a format for the storage of chunked, compressed, N-dimensional arrays inspired by HDF5, h5py and bcolz.
* [__OME-NGFF__](https://ngff.openmicroscopy.org) (Next Generation File Format) is a format based on zarr for the storage of biomedical imaging data.
* [__NIfTI__](https://nifti.nimh.nih.gov) (Neuroimaging Informatics Technology Initiative) is a single-file/single-resolution storage format for 3D+ neuroimaging data.
* [__BIDS__](https://bids-specification.readthedocs.io) (Brain Imaging Data Structure) is a simple and intuitive way to organize and describe data.


## Introduction

As biomedical imaging scales up, it is making more and more use of remote storage, remote computing and remote visualization.
Classical file formats—which store array data contiguously in a single file—are limited at large scales as
1. They often do not store data at multiple resolutions;
2. They do not offer efficient parallel access to data chunks.
   
These limits are very clear when it comes to visualizing very large data volumes, which cannot be loaded in full in memory.
In this context, it is preferable to only load the data required to display the current view (either a large field-of-view
at low-resolution, or a small field-of-view at high resolution).

The zarr format was developed to bypass the limitations of single-file formats such as HDF5. 
The microscopy community is currently developping its own standard for cloud-friendly biomedical imaging data: OME-NGFF. 
It build on zarr and adds rules for storing multi-resolutions images and medical-specific metadata such as axis names 
and voxel sizes. However, the microscopy community has needs in terms of metadata and coordinate-space description 
that are relatively complex, as they need to conform to different organs, different acquisition systems, or different tissue 
processing pipelines. This has drastically slowed down the adoption of a coordinate transform standard—the current version (0.5)
only handles canonical scales and offsets—and has also made the future coordinate transform standard much more complicated.

In contrast, the neuroimaging community has adopted and used a standard "world" coordinate frame for decades, where
```
+x = left      -> right
+y = posterior -> anterior
+z = inferior  -> superior
```
An affine transform is used to map from  the F-ordered voxel space (i, j, k) to world space (x, y, z). 
The neuroimaging community has also created a simple data exchange format—NIfTI—that has been widely embraced and 
is the mandatory file format in standardization efforts such as BIDS. However, the lack of multiresolution/chunk 
support in NIfTI has lead BIDS to adopt OME-TIFF and OME-ZARR as standard format for its microscopy component.

The NIfTI-Zarr (`nii.zarr`) specification attempts to merge the best of both worlds, in the simplest possible way. 
Like NIfTI, it aims to make the implementation of I/O libraries as simple as possible, to maximize chances that it 
gets adopted by the community. Its guiding principles are
* __OME-NGFF compliant:__ any `nii.zarr` file should be a valid `ome.zarr` file.
* __OME-NGFF minimal:__ only implements the minimum set of metadata necessary to describe [multi-resolution] neuroimaging data
* __NIfTI-complicant:__ the binary nifti header should be stored in the [group-level] `.zattrs` object.
* __NIfTI-priority:__ if metadata conflict across the nifti header and OME attributes, the nifti metadata should take precendence.

The simplicity of these guiding principles should make the adoption of `nii.zarr` in cloud environments (almost) as straightforward 
as the adoption of compressed-NIfTI (`.nii.gz`).


## Main differences with NIfTI and/or OME-NGFF

* Following the OME-NGFF specifcation, dimensions are ordered as [T, C, Z, Y, X] (in C order)
  as opposed to [C, T, Z, Y, X].
* To conform with the NIfTI expectation, on load data should be returned as a [C, T, Z, Y, X] array
  (in C order; [X, Y, Z, T, C] in F order).

### NIfTI features that are not supported by NIfTI-Zarr

* Any file with has more than 5 dimensions

### OME-Zarr features that are not supported by NIfTI-Zarr


### NIfTI > Zarr maps

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

In Zarr, the byte order **MUST** be specified by prepending one of `{"|", "<", ">"}` to the data type string.

## Format structure

### NIfTI header

The nifti header ([v1](https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h) or [v2](https://nifti.nimh.nih.gov/pub/dist/doc/nifti2.h)) 
should be encoded as a string using [base64](https://en.wikipedia.org/wiki/Base64) encoding and saved in the group-level `.zattrs` under 
the `"nifti"` key.

As a reminder, the nifti1 header has the following structure:

```C

                        /*************************/  /************************/
struct nifti_1_header { /* NIFTI-1 usage         */  /* ANALYZE 7.5 field(s) */
                        /*************************/  /************************/

                                           /*--- was header_key substruct ---*/
 int   sizeof_hdr;    /*!< MUST be 348           */  /* int sizeof_hdr;      */
 char  data_type[10]; /*!< ++UNUSED++            */  /* char data_type[10];  */
 char  db_name[18];   /*!< ++UNUSED++            */  /* char db_name[18];    */
 int   extents;       /*!< ++UNUSED++            */  /* int extents;         */
 short session_error; /*!< ++UNUSED++            */  /* short session_error; */
 char  regular;       /*!< ++UNUSED++            */  /* char regular;        */
 char  dim_info;      /*!< MRI slice ordering.   */  /* char hkey_un0;       */

                                      /*--- was image_dimension substruct ---*/
 short dim[8];        /*!< Data array dimensions.*/  /* short dim[8];        */
 float intent_p1 ;    /*!< 1st intent parameter. */  /* short unused8;       */
                                                     /* short unused9;       */
 float intent_p2 ;    /*!< 2nd intent parameter. */  /* short unused10;      */
                                                     /* short unused11;      */
 float intent_p3 ;    /*!< 3rd intent parameter. */  /* short unused12;      */
                                                     /* short unused13;      */
 short intent_code ;  /*!< NIFTI_INTENT_* code.  */  /* short unused14;      */
 short datatype;      /*!< Defines data type!    */  /* short datatype;      */
 short bitpix;        /*!< Number bits/voxel.    */  /* short bitpix;        */
 short slice_start;   /*!< First slice index.    */  /* short dim_un0;       */
 float pixdim[8];     /*!< Grid spacings.        */  /* float pixdim[8];     */
 float vox_offset;    /*!< Offset into .nii file */  /* float vox_offset;    */
 float scl_slope ;    /*!< Data scaling: slope.  */  /* float funused1;      */
 float scl_inter ;    /*!< Data scaling: offset. */  /* float funused2;      */
 short slice_end;     /*!< Last slice index.     */  /* float funused3;      */
 char  slice_code ;   /*!< Slice timing order.   */
 char  xyzt_units ;   /*!< Units of pixdim[1..4] */
 float cal_max;       /*!< Max display intensity */  /* float cal_max;       */
 float cal_min;       /*!< Min display intensity */  /* float cal_min;       */
 float slice_duration;/*!< Time for 1 slice.     */  /* float compressed;    */
 float toffset;       /*!< Time axis shift.      */  /* float verified;      */
 int   glmax;         /*!< ++UNUSED++            */  /* int glmax;           */
 int   glmin;         /*!< ++UNUSED++            */  /* int glmin;           */

                                         /*--- was data_history substruct ---*/
 char  descrip[80];   /*!< any text you like.    */  /* char descrip[80];    */
 char  aux_file[24];  /*!< auxiliary filename.   */  /* char aux_file[24];   */

 short qform_code ;   /*!< NIFTI_XFORM_* code.   */  /*-- all ANALYZE 7.5 ---*/
 short sform_code ;   /*!< NIFTI_XFORM_* code.   */  /*   fields below here  */
                                                     /*   are replaced       */
 float quatern_b ;    /*!< Quaternion b param.   */
 float quatern_c ;    /*!< Quaternion c param.   */
 float quatern_d ;    /*!< Quaternion d param.   */
 float qoffset_x ;    /*!< Quaternion x shift.   */
 float qoffset_y ;    /*!< Quaternion y shift.   */
 float qoffset_z ;    /*!< Quaternion z shift.   */

 float srow_x[4] ;    /*!< 1st row affine transform.   */
 float srow_y[4] ;    /*!< 2nd row affine transform.   */
 float srow_z[4] ;    /*!< 3rd row affine transform.   */

 char intent_name[16];/*!< 'name' or meaning of data.  */

 char magic[4] ;      /*!< MUST be "ni1\0" or "n+1\0". */

} ;                   /**** 348 bytes total ****/
```

### Single resolution

Directory structure:
```
└── mri.nii.zarr              # A nifti volume converted to Zarr.
    │
    ├── .zgroup               # An image is a Zarr array.
    ├── .zattrs               # Attributes are stored in the .zattrs file and include "nifti" (see below).
    │
    └── 0
        └─ t                  # Chunks are stored with the nested directory layout.
           └─ c               # All but the last chunk element are stored as directories.
              └─ z            # The terminal chunk is a file. Together the directory and file names
                 └─ y         # provide the "chunk coordinate" (t, c, z, y, x), where the maximum coordinate
                    └─ x      # will be dimension_size / chunk_size.

```

`.zattrs` file
```python
{
    # Nifti header encoded in base64
    "nifti": "...",
    # OME-NGFF resolution specification
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
            "datasets": [{
                 "path": "0",
                 "coordinateTransformations": [{
                     # the voxel size for the first scale level (0.5 millimeter)
                     "type": "scale",
                     "scale": [1.0, 1.0, 0.5, 0.5, 0.5]
                 }]
            }],
            "coordinateTransformations": [{
                # the time unit (1.0 seconds), which is the same for each scale level
                "type": "scale",
                "scale": [0.1, 1.0, 1.0, 1.0, 1.0]
            }],
        }
    ]
}
```

### Multi resolution

Directory structure:
```
└── mri.nii.zarr              # A nifti volume converted to Zarr.
    │
    ├── .zgroup               # Each image is a Zarr group, or a folder, of other groups and arrays.
    ├── .zattrs               # Group level attributes are stored in the .zattrs file and include
    │                         # "multiscales" and "nifti" (see below). In addition, the group level attributes
    │                         # may also contain "_ARRAY_DIMENSIONS" for compatibility with xarray if
    |                         # this group directly contains multi-scale arrays.
    │
    ├── 0                     # Each multiscale level is stored as a separate Zarr array,
    │   ...                   # which is a folder containing chunk files which compose the array.
    └── n                     # The name of the array is arbitrary with the ordering defined by
        │                     # by the "multiscales" metadata, but is often a sequence starting at 0.
        │
        ├── .zarray           # All image arrays must be up to 5-dimensional
        │                     # with the axis of type time before type channel, before spatial axes.
        │
        └─ t                  # Chunks are stored with the nested directory layout.
           └─ c               # All but the last chunk element are stored as directories.
              └─ z            # The terminal chunk is a file. Together the directory and file names
                 └─ y         # provide the "chunk coordinate" (t, c, z, y, x), where the maximum coordinate
                    └─ x      # will be dimension_size / chunk_size.

```

`.zattrs` file
```python
{
    # Nifti header encoded in base64
    "nifti": "...",
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
                        # the voxel size for the second scale level (downscaled by a factor of 2 -> 1 millimeter)
                        "type": "scale",
                        "scale": [1.0, 1.0, 1.0, 1.0, 1.0]
                    }]
                },
                {
                    "path": "2",
                    "coordinateTransformations": [{
                        # the voxel size for the third scale level (downscaled by a factor of 4 -> 2 millimeter)
                        "type": "scale",
                        "scale": [1.0, 1.0, 2.0, 2.0, 2.0]
                    }]
                }
            ],
            "coordinateTransformations": [{
                # the time unit (1.0 seconds), which is the same for each scale level
                "type": "scale",
                "scale": [0.1, 1.0, 1.0, 1.0, 1.0]
            }],
        }
    ]
}
```
