#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MicroCT image processing utilities.

"""

__author__ = 'Gianluca Iori'
__date_created__ = '2021-03-28'
__date__ = '2023-01-14'
__copyright__ = 'Copyright (c) 2022, ORMIR'
__docformat__ = 'restructuredtext en'
__license__ = "MIT"
__version__ = "1.4"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"

import numpy as np
import png
import os
import tifffile
import matplotlib.pyplot as plt

def touint8(data_3D, range=None, quantiles=None, numexpr=True):
    """Normalize and convert data to uint8.

    Parameters
    ----------
    data_3D
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

    def convert8bit():

        data_3D_float = data_3D.astype(np.float32, copy=False)
        df = np.float32(data_max - data_min)
        mn = np.float32(data_min)

        if numexpr:
            import numexpr as ne

            scl = ne.evaluate('0.5+255*(data_3D_float-mn)/df', truediv=True)
            ne.evaluate('where(scl<0,0,scl)', out=scl)
            ne.evaluate('where(scl>255,255,scl)', out=scl)
            return scl.astype(np.uint8)
        else:
            data_3D_float = 0.5 + 255 * (data_3D_float - mn) / df
            data_3D_float[data_3D_float < 0] = 0
            data_3D_float[data_3D_float > 255] = 255
            return np.uint8(data_3D_float)

    if range == None:

        # if quantiles is empty data is scaled based on its min and max values
        if quantiles == None:
            data_min = np.nanmin(data_3D)
            data_max = np.nanmax(data_3D)
            data_max = data_max - data_min
            return convert8bit()
        else:
            [data_min, data_max] = np.quantile(np.ravel(data_3D), quantiles)
            return convert8bit()

    else:
        # ignore quantiles input if given
        if quantiles is not None:
            print('quantiles input ignored.')

        data_min = range[0]
        data_max = range[1]
        return convert8bit()

def to01(data_3D):
    """Normalize data to 0-1 range.

    Parameters
    ----------
    data_3D
        Input data.

    Returns
    -------
    data_3D : float32
        Normalized data.
    """
    import numexpr as ne

    data_3D = data_3D.astype(np.float32, copy=False)
    data_min = np.nanmin(data_3D)
    data_max = np.nanmax(data_3D)
    df = np.float32(data_max - data_min)
    mn = np.float32(data_min)
    scl = ne.evaluate('(data_3D-mn)/df', truediv=True)
    return scl.astype(np.float32)

def writemidplanes(data_3D, fileout, slice_x=-1, slice_y=-1, slice_z=-1):
    """Plot orthogonal mid-planes through 3D dataset and save them as images.
    Uses pypng for writing .PNG files.

    Parameters
    ----------
    data
        Input 3D image data.
    fileout : str
        Output .PNG image file name.
    slice_x : int
        X-slice number.
    slice_y : int
        Y-slice number.
    slice_z : int
        Z-slice number.
    """

    if data_3D.ndim == 3:

        if slice_x == -1:
            slice_x = int(data_3D.shape[2] / 2)
        if slice_y == -1:
            slice_y = int(data_3D.shape[1] / 2)
        if slice_z == -1:
            slice_z = int(data_3D.shape[0] / 2)

        filename, ext = os.path.splitext(fileout)
        with open(filename + '_XY.png', 'wb') as midplaneXY:
            pngWriter = png.Writer(data_3D.shape[2], data_3D.shape[1], greyscale=True, alpha=False, bitdepth=8)
            pngWriter.write(midplaneXY, touint8(data_3D[int(slice_z), :, :]))

        with open(filename + '_XZ.png', 'wb') as midplaneXZ:
            pngWriter = png.Writer(data_3D.shape[2], data_3D.shape[0], greyscale=True, alpha=False, bitdepth=8)
            pngWriter.write(midplaneXZ, touint8(data_3D[:, int(slice_y), :]))

        with open(filename + '_YZ.png', 'wb') as midplaneYZ:
            pngWriter = png.Writer(data_3D.shape[1], data_3D.shape[0], greyscale=True, alpha=False, bitdepth=8)
            pngWriter.write(midplaneYZ, touint8(data_3D[:, :, int(slice_x)]))

def plot_midplanes(data_3D, slice_x=-1, slice_y=-1, slice_z=-1):
    """Plot orthogonal cross-sections through 3D dataset.

    Parameters
    ----------
    data_3D
        Input 3D image data.
    slice_x : int
        X-slice number.
    slice_y : int
        Y-slice number.
    slice_z : int
        Z-slice number.
    """

    if slice_x == -1:
        slice_x = int(data_3D.shape[2] / 2)
    if slice_y == -1:
        slice_y = int(data_3D.shape[1] / 2)
    if slice_z == -1:
        slice_z = int(data_3D.shape[0] / 2)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.imshow(data_3D[slice_z, :, :])
    ax2.imshow(data_3D[:, slice_y, :])
    ax3.imshow(data_3D[:, :, slice_x])

def plot_projections(data_3D, projection='max'):
    """Plot orthogonal projections of 3D dataset.

    Parameters
    ----------
    data_3D
        Input 3D image data.
    projection : str
        Projection method. Available choices are 'max', 'min'.
    """

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    if projection=='max':
        ax1.imshow(np.max(data_3D, 0))
        ax2.imshow(np.max(data_3D, 1))
        ax3.imshow(np.max(data_3D, 2))
    elif projection=='min':
        ax1.imshow(np.min(data_3D, 0))
        ax2.imshow(np.min(data_3D, 1))
        ax3.imshow(np.min(data_3D, 2))

def read_tiff_stack(filename, range=None, zfill=4):
    """Read stack of tiff files. Searches all files in parent folder and opens them as a stack of images.

    Parameters
    ----------
    filename
        One of the stack images.
    range : [int, int]
        Control load slices range.
    zfill : int
        Number of leading zeros in file names.

    TO DO:
    ----------
    - check that folder contains only .TIFF files; skip the rest
    """

    # search all files in parent folder; create filenames list
    stack_files = [os.path.join(os.path.dirname(filename), f) for f in os.listdir(os.path.dirname(filename))
                     if os.path.isfile(os.path.join(os.path.dirname(filename), f))]
    stack_files.sort()

    if range is not None:
        import re
        slice_in = [i for i, item in enumerate(stack_files) if re.search(str(range[0]).zfill(4) + ".", item)]
        slice_end = [i for i, item in enumerate(stack_files) if re.search(str(range[1]).zfill(4) + ".", item)]

        if len(slice_in) == 1 and len(slice_end) == 1:
            stack_files = stack_files[slice_in[0]:slice_end[0]]
        else:
            import warnings
            warnings.warn('Given slice range is ambiguous or non existing.. loading whole stack.')

    # load stack using tifffile
    return tifffile.imread(stack_files)

def bbox(bw, pad=0, dsize=None, verbose=None):
    """Bounding BOX limits of input binary image.

    Parameters
    ----------
    bw : bool
        Binary image.
    pad : int
        Add padding of given number of pixels to the BBOX limits.
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

        # add padding
        row0 = row0 - pad
        rowd = row1 - row0 + pad
        col0 = col0 - pad
        cold = col1 - col0 + pad
        slice0 = slice0 - pad
        sliced = slice1 - slice0 + pad

        if pad > 0:
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

def crop(data_3D, crop_origin, crop_size):
    """Crop 3D image given crop origin and size.

    Parameters
    ----------
    data_3D
        Input data.
    crop_origin : [int, int, int]
        Crop origin [Z,Y,X].
    crop_size : [int, int, int]
        Crop size [Z,Y,X].

    Returns
    -------
    output
        Cropped data.
    """
    return data_3D[crop_origin[2]:crop_origin[2] + crop_size[2], crop_origin[0]:crop_origin[0] + crop_size[0], crop_origin[1]:crop_origin[1] + crop_size[1]]