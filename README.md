# nifti-zarr draft specification

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

In contrast, the neuroimaging community has adopted and used a standard "world" coordinate frame for decades
([left->right, posterior->anterior, inferior->superior] and its variations). The neuroimaging community has 
also created a simple data exchange format--NIfTI--that has been widely embraced and is the mandatory file format
in standardization efforts such as BIDS.

The nifti-zarr (`nii.zarr`) specification attempts to merge the best of both worlds, in the simplest possible way. 
Like nifti, it aims to make the implementation of I/O libraries as simple as possible, to maximize chances that it 
gets adopted by the community. Its guiding principles are
* __Be OME-NGFF compliant:__ any `nii.zarr` file should be a valid `ome.zarr` file.
* __Be OME-NGFF minimal:__ only implements the minimum set of metadata necessary to describe [multi-resolution] neuroimaging data
* __Be nifti-complicant:__ the binary nifti header should be stored in the [group-level] `.zattrs` object.
* __Have nifti-priority:__ if metadata conflict across the nifti header and OME attributes, the nifti metadata should take precendence.
