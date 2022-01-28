#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Bone histomorphometry and image processing methods.

"""

__author__ = ['Gianluca Iori']
__date_created__ = '2021-11-03'
__date__ = '2021-11-03'
__copyright__ = 'Copyright (c) 2021, JC|MSK'
__docformat__ = 'restructuredtext en'
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"

import numpy as np
from scipy import ndimage
from skimage import measure, morphology
import logging
from tqdm import tqdm

def centerofmass(bwimage):
    """Center Of Mass (COM) of binary image.

    Parameters
    ----------
    bwimage: bool
        Binary image. Can be 2D and 3D.

    Returns
    -------
    cmassx_array
        X-coordinate array of the COM. If input is 3D, an array of the slicewise COMs is returned.
    cmassy_array
        Y-coordinate array of the COM.
    """

    if bwimage.ndim == 3:

        # output arrays initialization
        cmassx_array = np.zeros([bwimage.shape[0]])
        cmassy_array = np.zeros([bwimage.shape[0]])

        for slice in range(0, bwimage.shape[0]):
            y = np.sum(bwimage[slice,:,:], 1)
            cmassy = np.inner(y, np.arange(0, y.size))
            cmassy_array[slice] = cmassy / np.sum(y)

            x = np.sum(bwimage[slice, :, :], 0)
            cmassx = np.inner(x, np.arange(0, x.size))
            cmassx_array[slice] = cmassx / np.sum(x)

    elif bwimage.ndim == 2:
        y = np.sum(bwimage, 1)
        cmassy = np.inner(y, np.arange(0, y.size))
        cmassy_array = cmassy / np.sum(y)

        x = np.sum(bwimage, 0)
        cmassx = np.inner(x, np.arange(0, x.size))
        cmassx_array = cmassx / np.sum(x)

    return cmassx_array, cmassy_array

def remove_unconnected(bwimage):
    """Remove all unconnected voxels. Returns a binary of the largest connected cluster.

    Parameters
    ----------
    bwimage
        Binary image.

    Returns
    -------
    bwcluster
        Binary image of the largest connected cluster of voxels.
    """

    # label the BW image
    # [labels, n_labels] = measure.label(bwimage, None, True)
    [labels, n_labels] = measure.label(bwimage, None, True, 1)

    # count occurrences of each label
    occurrences = np.bincount(labels.reshape(labels.size))

    # find largest unconnected label
    largest_label_id = occurrences[1:].argmax() + 1
    bwcluster = labels == largest_label_id

    return bwcluster

def periosteummask(bwimage, closepixels=10, closevoxels=0, remove_objects_smaller_than=None, removeunconn=True, verbose=False):
    """Binary mask of periosteum (whole bone).

    Parameters
    ----------
    bwimage : bool
        Binary image. Can be 2D or 3D.
    closepixels : int
        Radius of DISK structuring element for 2D image closing.
    closevoxels : int
        Radius of CUBE structuring element for final 3D image closing.
    remove_objects_smaller_than : int
        Remove objects smaller than given size before periosteum mask calculation.
    removeunconn : bool
        Remove unconnected clusters of pixels/voxels from the calculated mask.
    verbose : bool
        Activate verbose output.

    Returns
    -------
    perimask : bool
        Binary mask of the whole bone (periosteum mask).
    """

    # verbose output
    if verbose:
        logging.basicConfig(level=logging.INFO)

    if remove_objects_smaller_than:
        logging.info('Preliminary removal of objects smaller than {} pixels.'.format(remove_objects_smaller_than))
        bwimage = morphology.remove_small_objects(bwimage, min_size=remove_objects_smaller_than)

    if bwimage.ndim == 3:

        # output arrays initialization
        perimask = np.zeros(np.shape(bwimage), dtype=bool)

        # 2D slice-wise imclose and fill
        logging.info('2D slice-wise image closing and filling.\n Structuring element DISK of radius: {}'.format(closepixels))

        for slice in tqdm(range(0, bwimage.shape[0])):
            perimask[slice,:,:] = ndimage.binary_fill_holes(morphology.binary_closing(bwimage[slice,:,:], morphology.disk(closepixels)))

        if removeunconn:
            # remove isolated clusters
            logging.info("Removing isolated clusters of voxels.")
            perimask = remove_unconnected(perimask)

        if closevoxels > 0:
            # final 3D imclose
            logging.info('Final 3D image closing.\n Structuring element CUBE of radius: {}'.format(closepixels))
            perimask = morphology.binary_closing(perimask, morphology.cube(closevoxels))

    elif bwimage.ndim == 2:
        # imclose and fill
        perimask = ndimage.binary_fill_holes(morphology.binary_closing(bwimage, morphology.disk(closepixels)))

        # remove isolated clusters
        if removeunconn:
            # remove isolated clusters
            logging.info("Removing isolated clusters of voxels.")
            perimask = remove_unconnected(perimask)

    return perimask