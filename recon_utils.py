#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Some reconstruction image processing utilities.

"""

"""
2DO:
- bounding box
- latex report

"""

__author__ = 'Gianluca Iori'
__date_created__ = '2021-03-28'
__date__ = '2021-08-27'
__copyright__ = 'Copyright (c) 2021, BEATS'
__docformat__ = 'restructuredtext en'
__license__ = "GPL"
__version__ = "0.0.2"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"

import numpy as np
import png
import os
import dxchange
import tifffile
# import matplotlib
import matplotlib.pyplot as plt

def touint8(data, quantiles=None):
    # scale data to uint8
    # if quantiles is empty data is scaled based on its min and max values
    if quantiles == None:
        data_min = np.min(data)
        data_max = np.max(data)
        data_max = data_max - data_min
        data = 255 * ((data - data_min) / data_max)
        return np.uint8(data)
    else:
        [q0, q1] = np.quantile(np.ravel(data), quantiles)
        q1 = q1 - q0
        data = 255 * ((data - q0) / q1)
        return np.uint8(data)


def writemidplanes(data, filename_out):
    # write orthogonal mid-planes through 3D dataset
    if data.ndim == 3:
        filename, ext = os.path.splitext(filename_out)
        with open(filename + '_XY.png', 'wb') as midplaneXY:
            pngWriter = png.Writer(data.shape[2], data.shape[1], greyscale=True, alpha=False, bitdepth=8)
            pngWriter.write(midplaneXY, touint8(data[int(data.shape[0] / 2), :, :]))

        with open(filename + '_XZ.png', 'wb') as midplaneXZ:
            pngWriter = png.Writer(data.shape[2], data.shape[0], greyscale=True, alpha=False, bitdepth=8)
            pngWriter.write(midplaneXZ, touint8(data[:, int(data.shape[1] / 2), :]))

        with open(filename + '_YZ.png', 'wb') as midplaneYZ:
            pngWriter = png.Writer(data.shape[1], data.shape[0], greyscale=True, alpha=False, bitdepth=8)
            pngWriter.write(midplaneYZ, touint8(data[:, :, int(data.shape[2] / 2)]))

def writemidplanesDxchange(data, filename_out):
    if data.ndim == 3:
        filename, ext = os.path.splitext(filename_out)
        dxchange.writer.write_tiff(touint8(data[int(data.shape[0] / 2), :, :]), fname=filename+'_XY.tiff', dtype='uint8')
        dxchange.writer.write_tiff(touint8(data[:, int(data.shape[1] / 2), :]), fname=filename + '_XZ.tiff', dtype='uint8')
        dxchange.writer.write_tiff(touint8(data[:, :, int(data.shape[2] / 2)]), fname=filename + '_YZ.tiff', dype='uint8')

def plot_midplanes(data_3D, slice_x=-1, slice_y=-1, slice_z=-1):
    # Plot midplanes of 3D data
    if slice_x == -1:
        slice_x = int(data_3D.shape[1] / 2)
    if slice_y == -1:
        slice_y = int(data_3D.shape[2] / 2)
    if slice_z == -1:
        slice_z = int(data_3D.shape[0] / 2)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.imshow(data_3D[slice_z, :, :])
    ax2.imshow(data_3D[:, slice_x, :])
    ax3.imshow(data_3D[:, :, slice_y])

def read_tiff_stack(filename):
    # Read a stack of tiffs from single slice filename.
    # Searches all files in parent folder and opens them as a stack of images.
    # TO DO:
    #     - check that folder contains only .TIFF files; skip the rest

    # search all files in parent folder; create filenames list
    tifffiles = [os.path.join(os.path.dirname(filename), f) for f in os.listdir(os.path.dirname(filename))
                     if os.path.isfile(os.path.join(os.path.dirname(filename), f))]
    tifffiles.sort()

    # load stack using tifffile
    return tifffile.imread(tifffiles)

def add_cap(data_3D, cap_thickness, cap_val):
    # Add caps of voxels with given GV to the input 3D data.
    # Caps are added on both ends along the Z-direction (first dataset dimension).
    data_3D_cap = np.ones([data_3D.shape[0]+2*cap_thickness, data_3D.shape[1], data_3D.shape[2]], data_3D.dtype)*cap_val
    data_3D_cap[cap_thickness:-cap_thickness, :, :] = data_3D
    return data_3D_cap
