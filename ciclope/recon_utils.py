#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tomographic reconstruction image processing utilities.

"""

"""
2DO:
- 2D bbox
- fast touint8
- latex report

"""

__author__ = 'Gianluca Iori'
__date_created__ = '2021-03-28'
__date__ = '2021-11-09'
__copyright__ = 'Copyright (c) 2021, JC|MSK'
__docformat__ = 'restructuredtext en'
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"

import numpy as np
import numexpr as ne
import png
import os
import dxchange
import tifffile
import matplotlib.pyplot as plt

def convert8bit(rec, data_min, data_max, numexpr=True):
    rec = rec.astype(np.float32, copy=False)
    df = np.float32(data_max-data_min)
    mn = np.float32(data_min)

    if numexpr:
        scl = ne.evaluate('0.5+255*(rec-mn)/df', truediv=True)
        ne.evaluate('where(scl<0,0,scl)', out=scl)
        ne.evaluate('where(scl>255,255,scl)', out=scl)
        return scl.astype(np.uint8)
    else:
        rec = 0.5+255*(rec-mn)/df
        rec[rec<0]=0
        rec[rec>255]=255
        return np.uint8(rec)

def touint8(data, range=None, quantiles=None, numexpr=True):
    """Normalize and convert data to uint8.

        Parameters
        ----------
        data
            Input data.
        range : [float, float]
            Control range for data normalization.
        quantiles : [float, float]
            Define range for data normalization through input data quantiles. If range is given this input is ignored.
        numexpr : bool
            Use fast numerical expression evaluator for NumPy (memory expensive).

        Returns
        -------
        output : uint8
            Normalized data.
        """

    if range == None:

        # if quantiles is empty data is scaled based on its min and max values
        if quantiles == None:
            data_min = np.nanmin(data)
            data_max = np.nanmax(data)
            data_max = data_max - data_min
            return convert8bit(data, data_min, data_max, numexpr)
        else:
            [q0, q1] = np.quantile(np.ravel(data), quantiles)
            return convert8bit(data, q0, q1, numexpr)

    else:
        # ignore quantiles input if given
        if quantiles is not None:
            print('quantiles input ignored.')

        return convert8bit(data, range[0], range[1], numexpr)

def to01(I):
    """Normalize data to 0-1 range.

    Parameters
    ----------
    I
        Input data.

    Returns
    -------
    I : float32
        Normalized data.
    """

    I = I.astype(np.float32, copy=False)
    data_min = np.nanmin(I)
    data_max = np.nanmax(I)
    df = np.float32(data_max - data_min)
    mn = np.float32(data_min)
    scl = ne.evaluate('(I-mn)/df', truediv=True)
    return scl.astype(np.float32)

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

def bbox(bw, offset=0, dsize=None, verbose=None):
    """Bounding BOX limits of input binary image.

    Parameters
    ----------
    bw : bool
        Binary image.
    offset : int
        Add offset of given number of pixels to the BBOX limits.
    dsize : int
        perform image close with disk structuring element of radius 'dsize' before calculating the BBOX.
    verbose
        Activate verbose graphical output

    Returns
    -------
    bbox_origin: int
        Origin [row col (slice)] of the BBOX inscribing True values in input image bw.
    bbox_size: int
        BBOX size [s_row s_col (s_slice)].
    """

    # DSIZE: remove artefacts > erode/dilate
    if dsize:
        raise IOError('dsize method not implemented yet.')

    if bw.ndim == 3:
        # project along each dimension
        maxROW = np.max(np.max(bw, 0), 1)
        maxCOL = np.max(np.max(bw, 0), 0)
        maxSLICE = np.max(np.max(bw, 1), 1)

        # find first and last True occurrences
        row0 = list(maxROW).index(True)
        row1 = len(maxROW) - list(maxROW[::-1]).index(True) - 1

        col0 = list(maxCOL).index(True)
        col1 = len(maxCOL) - list(maxCOL[::-1]).index(True) - 1

        slice0 = list(maxSLICE).index(True)
        slice1 = len(maxSLICE) - list(maxSLICE[::-1]).index(True) - 1

        # add offset
        row0 = row0 - offset
        rowd = row1 - row0 + offset
        col0 = col0 - offset
        cold = col1 - col0 + offset
        slice0 = slice0 - offset
        sliced = slice1 - slice0 + offset

        if offset > 0:
            # check if bbox exceeds image size
            if row0 < 0:
                row0=0
            if col0 < 0:
                col0=0
            if slice0 < 0:
                slice0=0

            bw_size = bw.shape
            if slice0 + sliced > bw_size[0]:
                sliced = bw_size[0] - slice0
            if row0 + rowd > bw_size[1]:
                rowd = bw_size[1] - row0
            if col0 + cold > bw_size[2]:
                cold = bw_size[2]-col0

        if verbose:
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.imshow(np.max(bw, 0))
            ax1.plot([col0, col0], [0, bw.shape[1]-1], 'r')
            ax1.plot([col0+cold, col0+cold], [0, bw.shape[1]-1], 'r')
            ax1.plot([0, bw.shape[2] - 1], [row0, row0], 'r')
            ax1.plot([0, bw.shape[2] - 1], [row0 + rowd, row0 + rowd], 'r')

            ax2.imshow(np.max(bw, 1))
            ax2.plot([col0, col0], [0, bw.shape[1] - 1], 'r')
            ax2.plot([col0 + cold, col0 + cold], [0, bw.shape[1] - 1], 'r')
            ax2.plot([0, bw.shape[0] - 1], [slice0, slice0], 'r')
            ax2.plot([0, bw.shape[0] - 1], [slice0 + sliced, slice0 + sliced], 'r')

        bbox_origin = [row0, col0, slice0]
        bbox_size = [rowd, cold, sliced]

        return bbox_origin, bbox_size

    if bw.ndim == 2:
        raise IOError('bbox method for 2D images not implemented yet.')

def crop(I, crop_origin, crop_size):
    return I[crop_origin[2]:crop_origin[2] + crop_size[2], crop_origin[0]:crop_origin[0] + crop_size[0], crop_origin[1]:crop_origin[1] + crop_size[1]]
