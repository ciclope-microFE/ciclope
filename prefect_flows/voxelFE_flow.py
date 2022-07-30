#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for prefect flow generation
"""

import numpy as np
from ciclope.utils.preprocess import resample, segment, remove_unconnected
from ciclope.core.voxelFE import vol2voxelfe, matpropdictionary
import ciclope.utils.recon_utils as ru
from skimage.filters import gaussian
import prefect
from prefect import task, Flow, Parameter

@task
def load_tiff_stack_prefect(filein: str) -> np.ndarray:
    logger = prefect.context.get("logger")
    logger.info("Reading image: %s" % filein)
    return ru.read_tiff_stack(filein)

@task
def smooth_prefect(image: np.ndarray, smooth: float) -> np.ndarray:
    if smooth != 0:
        logger = prefect.context.get("logger")
        logger.info("Apply gaussian smooth")
        image = gaussian(image, sigma=smooth, preserve_range=True)

    return image

@task
def resample_prefect(image: np.ndarray, voxelsize, resampling_factor=1):
    if (resampling_factor != 0) & (resampling_factor != 1):
        logger = prefect.context.get("logger")
        logger.info("Resampling input dataset with factor 1/{}".format(resampling_factor))
        image, voxelsize = resample(image, voxelsize, resampling_factor)
        logger.info("New voxelsize: {}".format(voxelsize))

    return image

@task
def update_voxelsize(voxelsize: float, resampling_factor: float) -> float:
    if (resampling_factor != 0) & (resampling_factor != 1):
        voxelsize*=resampling_factor
    return voxelsize

@task
def add_cap_prefect(image: np.ndarray, cap_thickness: int, cap_val: int) -> np.ndarray:
    if cap_thickness != 0:
        # add cap so that the mesh is closed
        logger = prefect.context.get("logger")
        logger.info("Adding cap to image")
        image = ru.add_cap(image, cap_thickness=cap_thickness, cap_val=cap_val)

    return image

@task
def segment_prefect(image: np.ndarray, threshold_value: float) -> np.ndarray:
    logger = prefect.context.get("logger")
    logger.info("Segmenting the dataset...")
    if threshold_value is None:
        logger.info("Use Otsu's method.")
    BW, T = segment(image, threshold_value)
    return BW

@task
def remove_unconnected_prefect(bwimage: np.ndarray) -> np.ndarray:
    logger = prefect.context.get("logger")
    logger.info("Removing unconnected clusters of voxels...")
    return remove_unconnected(bwimage)

@task
def generate_voxelfe_prefect(bwimage: np.ndarray, templatefile: str, filenameout: str, voxelsize: list):
    """voxelFE of binary volume data; the material property definition is assumed to be in the analysis template file
    """
    logger = prefect.context.get("logger")
    logger.info("Generating voxelFE model file with constant material property.")
    vol2voxelfe(bwimage, templatefile, filenameout, keywords=['NSET', 'ELSET'], voxelsize=voxelsize)

@task
def generate_voxelfe_mapping_prefect(image: np.ndarray, bwimage: np.ndarray, templatefile: str, filenameout: str, voxelsize: list, mapping: list):
    """
    # voxelFE of greyscale volume data; material mapping on

    :param bwimage:
    :param bwimage:
    :param templatefile:
    :param filenameout:
    :param voxelsize:
    :return:
    """
    # build dictionary of material properties
    matprop = matpropdictionary(mapping)

    # mask the data (keep only the largest connected volume)
    image[~bwimage] = 0

    vol2voxelfe(image, templatefile, filenameout, matprop, keywords=['NSET', 'ELSET', 'PROPERTY'], voxelsize=voxelsize)

def ciclope_voxelfe_flow(master_file: str):
    """**ciclope** prefect flow for voxelFE model generation.
    Maps and loops through a given list of voxelFE ciclopes.
    :return:
    """
    import pandas as pd
    with Flow("ciclope-voxelFE-flow") as flow:

        # input parameters
        filename = Parameter('filein', required=True)
        fileout = Parameter('fileout', required=True)
        s = Parameter('smooth')
        vs = Parameter('vs')
        r = Parameter('r')
        t = Parameter('t')
        cap_t = Parameter('cap_t')
        cap_val = Parameter('cap_val')
        template = Parameter('template')

        # ciclope flow
        I = load_tiff_stack_prefect.map(filename)
        I = smooth_prefect.map(image=I, smooth=s)
        I = resample_prefect.map(image=I, voxelsize=vs, resampling_factor=r)
        vs = update_voxelsize.map(voxelsize=vs, resampling_factor=r)
        I = add_cap_prefect.map(image=I, cap_thickness=cap_t, cap_val=cap_val)
        BW = segment_prefect.map(image=I, threshold_value=t)
        BW = remove_unconnected_prefect.map(bwimage=BW)
        generate_voxelfe_prefect.map(bwimage=BW, templatefile=template, filenameout=fileout, voxelsize=vs)

    # load master table
    df = pd.read_csv(master_file)

    # set all empty cells to None
    df = df.where(pd.notnull(df), None)

    # flow = ciclope_flow()

    parameters = df[df['run'] == 1].to_dict('list')
    del parameters['run']
    # parameters2 = {'filein': df['filein'][df['run'] == 1].tolist(),
    #                'smooth': df['smooth'][df['run'] == 1].tolist(),
    #                'vs': df['vs'][df['run'] == 1].tolist(),
    #                'r': df['r'][df['run'] == 1].tolist(),
    #                't': df['t'][df['run'] == 1].tolist(),
    #                'fileout': df['fileout'][df['run'] == 1].tolist(),
    #                'cap_t': df['cap_t'][df['run'] == 1].tolist(),
    #                'cap_val': df['cap_val'][df['run'] == 1].tolist(),
    #                'template': df['template'][df['run'] == 1].tolist(),
    #                }

    state = flow.run(parameters)

    type(state._result.value)
    task_ref = flow.get_tasks()[1]
    ru.plot_midplanes(state.result[task_ref]._result.value[0])

    return flow
