#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ciclope - Image preprocessing Module
"""

from scipy import ndimage, misc
from skimage.filters import threshold_otsu, gaussian

def segment(image, threshold_value):
    """Threshold image.

    Parameters
    ----------
    image
        Image data.
    threshold_value (optional)
        Threshold value. If empty an Otsu threshold is calculated.

    Returns
    -------
    BWimage
        Binary image after thresholding.
    T
        Threshold value.
    """

    # we do want bone = 1 and background = 0;
    if threshold_value is None:
        # use Otsu if threshold input not specified
        T = threshold_otsu(image)

    else:
        T = int(threshold_value)

    # apply the threshold
    return image > T, T

def resample(image, voxelsize, resampling_factor):
    """Resize image.

    Parameters
    ----------
    image
        Image data.
    voxelsize
        Voxel size.
    resampling_factor
        Scaling factor.

    Returns
    -------
    image
        Resized image.
    voxelsize
        Voxel size after rescaling.
    """

    # resize the 3D data using spline interpolation of order 2
    image = ndimage.zoom(image, 1 / resampling_factor, output=None, order=2)

    # correct voxelsize
    voxelsize = voxelsize * resampling_factor

    return image, voxelsize
