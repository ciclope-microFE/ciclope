#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Image preprocessing Module
"""

from scipy import ndimage, misc
from skimage.filters import threshold_otsu, gaussian

def segment(image, threshold_value):
    # we do want bone = 1 and background = 0;
    if threshold_value is None:
        # use Otsu if threshold input not specified
        T = threshold_otsu(image)

    else:
        T = int(threshold_value)

    # apply the threshold
    return image > T, T

def resample(image, voxelsize, resampling_factor):
    # resize the 3D data using spline interpolation of order 2
    image = ndimage.zoom(image, 1 / resampling_factor, output=None, order=2)

    # correct voxelsize
    voxelsize = voxelsize * resampling_factor

    return image, voxelsize

