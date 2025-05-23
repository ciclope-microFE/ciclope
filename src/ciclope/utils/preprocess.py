#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ciclope image preprocessing module
"""

import numpy as np
from scipy import ndimage
from skimage import measure, morphology
from skimage.filters import threshold_otsu
import logging
from tqdm import tqdm
from . import recon_utils as ru
import cv2
import os
import imageio
from PIL import Image
import math

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

    return labels == largest_label_id

def fill_voids(I, fill_val=None, makecopy=False):
    """Fill voids within color image with given value.

    Parameters
    ----------
    I
        Input color image.
    fill_val
        Filling value.
    makecopy : bool
        Make copy of input image.

    Returns
    -------
    I_filled
        Filled image.
    """

    if fill_val is None:
        fill_val = I.max()

    # binarize and label inverse of the input image
    [labels, n_labels] = measure.label(~(I>0), None, True, 1)

    # count occurrences of each label
    occurrences = np.bincount(labels.reshape(labels.size))

    # find and delete largest label (background)
    largest_label_id = occurrences[1:].argmax() + 1
    labels[labels == largest_label_id] = 0

    if makecopy:
        I_filled = I.copy()
        I_filled[labels != 0] = fill_val
        return I_filled
    else:
        I[labels != 0] = fill_val
        return I

def remove_largest(bwimage):
    """Remove largest cluster of voxels in binary image.

    Parameters
    ----------
    bwimage
        Binary image.

    Returns
    -------
    bwcluster
        Binary image in which the largest cluster of voxels is removed.
    """

    # label the BW image
    [labels, n_labels] = measure.label(bwimage, None, True, 1)

    # count occurrences of each label
    occurrences = np.bincount(labels.reshape(labels.size))

    # find largest unconnected label
    largest_label_id = occurrences[1:].argmax() + 1
    labels[labels == largest_label_id] = 0

    return labels != 0

def add_cap(I, cap_thickness, cap_val):
    """Add caps to 3D image.
    Caps are added on both ends along the Z-direction (first dataset dimension). The thickness and color (Grey Value) of the added caps can be specified.

    Parameters
    ----------
    I
        3D data. Zeroes as background.
    cap_thickness : int
        Cap thickness in pixels.
    cap_val : float
        Cap grey value.

    Returns
    ----------
    I_cap
        Image with caps added.
    """

    I_cap = np.ones([I.shape[0]+2*cap_thickness, I.shape[1], I.shape[2]], I.dtype)*cap_val
    I_cap[cap_thickness:-cap_thickness, :, :] = I
    return I_cap

def embed(I, embed_depth, embed_dir, embed_val=None, pad=0, makecopy=False):
    """Add embedding to 3D image.
    Direction and depth of the embedded region should be given. Zeroes in the input image is considered to be background.

    Parameters
    ----------
    I
        3D data. Zeroes as background.
    embed_depth : int
        Embedding depth in pixels.
    embed_dir : str
        Embedding direction. Can be "-x", "+x", "-y", "+y", "-z", or "+z".
    embed_val : float
        Embedding grey value.
    pad = int
        Padding around bounding box of embedded area.
    makecopy : bool
        Make copy of the input image.

    Returns
    ----------
    I
        Embedded image. Same size as the input one.
    BW_embedding
        BW mask of the embedding area.
    """

    if embed_val is None:
        embed_val = I.max() + 1

    # binarize the input image
    BW_I = np.zeros(I.shape, dtype='bool')
    BW_I[I>0] = True

    # init embedding mask
    BW_embedding = np.zeros(BW_I.shape, dtype='bool')

    if embed_dir == "-z":
        dir = -1

        # start the embedding at first non-zero voxel
        embed_start = np.where(np.max(BW_I.max(1), 1) == True)[-1][-1]

        # project embedded area and find size of embedding
        bbox_origin, bbox_size = ru.bbox(BW_I[embed_start + (dir * embed_depth):, :, :], pad)

        # create embedding mask
        BW_embedding[embed_start + (dir * embed_depth):, bbox_origin[0]:bbox_origin[0] + bbox_size[0], bbox_origin[1]:bbox_origin[1] + bbox_size[1]] = True

    elif embed_dir == "+z":
        dir = 1

        # start the embedding at first non-zero voxel
        embed_start = np.where(np.max(BW_I.max(1), 1) == True)[0][0]

        # project embedded area and find size of embedding
        bbox_origin, bbox_size = ru.bbox(BW_I[:embed_start + (dir * embed_depth), :, :], pad)

        # create embedding mask
        BW_embedding[:embed_start + (dir * embed_depth), bbox_origin[0]:bbox_origin[0] + bbox_size[0], bbox_origin[1]:bbox_origin[1] + bbox_size[1]] = True

    elif embed_dir == "-x":
        dir = -1

        # start the embedding at first non-zero voxel
        embed_start = np.where(np.max(BW_I.max(0), 0) == True)[-1][-1]

        # project embedded area and find size of embedding
        bbox_origin, bbox_size = ru.bbox(BW_I[:, :, embed_start + (dir * embed_depth):], pad)

        # create embedding mask
        BW_embedding[bbox_origin[2]:bbox_origin[2] + bbox_size[2], bbox_origin[0]:bbox_origin[0] + bbox_size[0], embed_start + (dir * embed_depth):] = True

    elif embed_dir == "+x":
        dir = +1

        # start the embedding at first non-zero voxel
        embed_start = np.where(np.max(BW_I.max(0), 0) == True)[0][0]

        # project embedded area and find size of embedding
        bbox_origin, bbox_size = ru.bbox(BW_I[:, :, :embed_start + (dir * embed_depth)], pad)

        # create embedding mask
        BW_embedding[bbox_origin[2]:bbox_origin[2] + bbox_size[2], bbox_origin[0]:bbox_origin[0] + bbox_size[0], :embed_start + (dir * embed_depth)] = True

    elif embed_dir == "+y":
        dir = -1

        # start the embedding at first non-zero voxel
        embed_start = np.where(np.max(BW_I.max(0), 1) == True)[-1][-1]

        # project embedded area and find size of embedding
        bbox_origin, bbox_size = ru.bbox(BW_I[:, embed_start + (dir * embed_depth):, :], pad)

        # create embedding mask
        BW_embedding[bbox_origin[2]:bbox_origin[2] + bbox_size[2], embed_start + (dir * embed_depth):, bbox_origin[1]:bbox_origin[1] + bbox_size[1]] = True

    elif embed_dir == "-y":
        dir = +1

        # start the embedding at first non-zero voxel
        embed_start = np.where(np.max(BW_I.max(0), 1) == True)[0][0]

        # project embedded area and find size of embedding
        bbox_origin, bbox_size = ru.bbox(BW_I[:, :embed_start + (dir * embed_depth), :])

        # create embedding mask
        BW_embedding[bbox_origin[2]:bbox_origin[2] + bbox_size[2], :embed_start + (dir * embed_depth), bbox_origin[1]:bbox_origin[1] + bbox_size[1]] = True

    else:
        raise IOError("EMBED_DIR parameter unknown. Valid entries are -x, +x, -y, +y, -z, and +z.")

    # emboss embedding mask with the masked input image
    BW_embedding = remove_unconnected(BW_embedding & ~BW_I)

    # assign embedding val to input image
    if makecopy:
        I_output = I.copy()
        I_output[BW_embedding] = embed_val
        return I_output, BW_embedding
    else:
        I[BW_embedding] = embed_val
        return I, BW_embedding

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

def _collect_from_elements(boundEls, nodesSet, dirs):
    """
    Expands voxel elements into their 8 corner nodes and adds direction-tagged nodes to a set.

    Each element in `boundEls` is treated as a voxel defined by its lower corner `(z, y, x)`.
    This function computes all 8 corner node positions for each voxel and adds them to `nodesSet`
    with each specified direction from `dirs`.

    Parameters
    ----------
    boundEls : iterable of tuple of int
        An iterable of 3D indices representing the lower-front-left corner of each voxel element.
        Each element should be a tuple of three integers (z, y, x).

    nodesSet : set of tuple
        A set to which the resulting nodes will be added. Each node is represented as a 4-tuple
        `(z, y, x, dir)`, where `(z, y, x)` is the node position and `dir` is a directional tag.

    dirs : iterable
        A list or iterable of direction tags to associate with each corner node. Can be any
        hashable type (e.g., int, str).

    Returns
    -------
    None
        This function modifies `nodesSet` in-place and does not return a value.
    """
    for el in boundEls:
        z, y, x = el  # element index triple
        for dz in (0, 1):
            for dy in (0, 1):
                for dx in (0, 1):
                    # Add a tuple (z+dz, y+dy, x+dx, dir)
                    for d in dirs:
                        nodesSet.add((z + dz, y + dy, x + dx, d))

"""
Ciclope functions for the modelling of trabecular bone sample with a cylindrical shape
"""

## Cropping and Loading input data

def invert_images(input_folder, output_folder):
    """
    Invert pixel values in images and save them.

    Parameters:
    - input_folder (str): The path to the input folder containing the images.
    - output_folder (str): The path to the output folder where inverted images will be saved.

    Returns:
    None
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".bmp"):
            input_path = os.path.join(input_folder, filename)
            img = cv2.imread(input_path, cv2.IMREAD_GRAYSCALE)

            # Perform pixel value inversion
            inverted_img = cv2.bitwise_not(img)

            output_path = os.path.join(output_folder, filename)
            cv2.imwrite(output_path, inverted_img)
            

def convert_bmp_to_tiff(input_folder, tiff_folder):
    """
    Convert BMP images to TIFF format and save them.

    Parameters:
    - input_folder (str): The path to the input folder containing BMP images.
    - tiff_folder (str): The path to the output folder where TIFF images will be saved.

    Returns:
    - str: The path to the folder where TIFF images are saved.
    """
    if not os.path.exists(tiff_folder):
        os.makedirs(tiff_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".bmp"):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(tiff_folder, os.path.splitext(filename)[0] + ".tif")
            img = cv2.imread(input_path, cv2.IMREAD_GRAYSCALE)
            imageio.imsave(output_path, img)

    return tiff_folder


def crop_and_resize_images(input_folder, output_folder, diameter, pixel_spacing_mm=0.0195, res_X=257, res_Y=257):
    """
    Crop and resize images based on a provided diameter.

    Parameters:
    - input_folder (str): The path to the input folder containing the images.
    - output_folder (str): The path to the output folder where cropped and resized images will be saved.
    - diameter (float): The diameter in millimeters for cropping and resizing.
    - pixel_spacing_mm (float, optional): The pixel spacing in millimeters, used to convert the diameter into pixels.
    - res_X (int, optional): The target width in pixels for resizing.
    - res_Y (int, optional): The target height in pixels for resizing.

    Returns:
    None
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Calculate the radius in pixels based on the provided diameter
    radius_in_pixels = int(diameter / (2 * pixel_spacing_mm))

    for filename in os.listdir(input_folder):
        if filename.endswith((".tif", ".bmp")):
            input_path = os.path.join(input_folder, filename)
            img = cv2.imread(input_path, cv2.IMREAD_GRAYSCALE)

            h, w = img.shape

            y, x = np.ogrid[:h, :w]
            mask = (x - w//2)**2 + (y - h//2)**2 <= radius_in_pixels**2

            cropped_img = img.copy()
            cropped_img[~mask] = 0

            non_zero_x = np.where(cropped_img.any(axis=0))
            non_zero_y = np.where(cropped_img.any(axis=1))
            left, right = non_zero_x[0][0], non_zero_x[0][-1]
            top, bottom = non_zero_y[0][0], non_zero_y[0][-1]

            cropped_img = cropped_img[top:bottom+1, left:right+1]
            resized_img = cv2.resize(cropped_img, (res_X, res_Y), interpolation=cv2.INTER_LINEAR)

            output_path = os.path.join(output_folder, filename)
            Image.fromarray(resized_img).save(output_path)
            
def crop_images(input_folder, output_folder, diameter, pixel_spacing_mm=0.0195):
    """
    Crop images based on a provided diameter.

    Parameters:
    - input_folder (str): The path to the input folder containing the images.
    - output_folder (str): The path to the output folder where cropped images will be saved.
    - diameter (float): The diameter in millimeters for cropping.
    - pixel_spacing_mm (float, optional): The pixel spacing in millimeters, used to convert the diameter into pixels. Default is 0.0195 mm.

    Returns:
    None
    """
    if not os.path.exists(output_folder):
        print(f"Creating output folder: {output_folder}")
        os.makedirs(output_folder)
    
    # Calcolo del raggio in pixel basato sul diametro fornito
    radius_in_pixels = int(diameter / (2 * pixel_spacing_mm))

    for filename in os.listdir(input_folder):
        if filename.endswith((".tif", ".bmp")):
            input_path = os.path.join(input_folder, filename)
            img = cv2.imread(input_path, cv2.IMREAD_GRAYSCALE)

            h, w = img.shape

            # Creazione della maschera circolare
            y, x = np.ogrid[:h, :w]
            mask = (x - w//2)**2 + (y - h//2)**2 <= radius_in_pixels**2

            # Applichiamo la maschera all'immagine
            cropped_img = img.copy()
            cropped_img[~mask] = 0

            # Salvataggio dell'immagine risultante
            output_path = os.path.join(output_folder, filename)
            cv2.imwrite(output_path, cropped_img)
            
def replace_ElType_ref(filename, old_word, new_word):
    """
    Replace element type word in a `.inp` file for correct input to Calculix.

    This function reads a `.inp` file, replaces all instances of the specified
    old element type word with a new word, and writes the updated content
    back to the same file. This is useful for modifying finite element models
    for compatibility with Calculix or other finite element software.

    Parameters
    ----------
    filename : str
        The path to the `.inp` file that needs modification. This should include
        the full file name and path.
    old_word : str
        The word (typically an element type identifier) that is to be replaced 
        in the file.
    new_word : str
        The word to replace the old word with. This should be the new element
        type or identifier that is compatible with Calculix or the desired software.

    Returns
    -------
    None
    """

    # Read the content of the file
    with open(filename, 'r') as file:
        file_content = file.read()

    # Replace the old word with the new word
    new_content = file_content.replace(old_word, new_word)

    # Write the modified content back to the file
    with open(filename, 'w') as file:
        file.write(new_content)