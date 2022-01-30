#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for voxel Finite Element model generation
"""

import os
import logging
import numpy as np
from datetime import datetime

def vol2voxelfe(voldata, templatefile, fileout, matprop=None, keywords=['NSET', 'ELSET'], voxelsize=[1, 1, 1], eltype='C3D8', matpropbits=8, refnode=None, verbose=False):
    """Generate ABAQUS voxel Finite Element (FE) input file from 3D volume data.
    The file written is an input file (.INP) in ABAQUS syntax that can be solved using ABAQUS or CALCULIX.
    The user can define a material mapping strategy for the conversion of local GVs to local material properties in the FE model.
    Material mapping laws are defined in separate template file(s) (see "prop.inp" and "property_temp_bone.inp" for examples).
    Boundary conditions, analysis type and output requests are defined in a separate template file (see "tmp.inp" for an example).
    Info on analysis definition at: https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-m-Sim-sb.htm#simacae-m-Sim-sb

    Parameters
    ----------
    voldata : ndarray
        3D voxel data.
    templatefile : str
        Analysis template file.
    fileout : str
        Output .INP file.
    matprop : dict
        Dictionary of material properties for material property mapping:
        matprop = {
            "file": ["prop.inp", "property_temp_bone.inp", ...],
            "range": [[250, 255], [0, 250], ...],
        }
    keywords : str
        SUPPORTED ABAQUS KEYWORDS:

        For a list of all Abaqus keywords and their description visit:
        https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-c-gen-kwbrowser.htm#simacae-c-gen-kwbrowser__simacae-gen-xsl-U

        * 'NODE':
                Specify nodal coordinates. (This cannot be turned off)
        * 'ELEMENT':
                Define elements by giving their nodes. (This cannot be turned off)
        * 'NSET':
                Create boundary node sets. (Default = ON)
                If 'NSET' is specified, the following node sets are created:
                  - NODES_Y0: Nodes on SOUTH (Y-) surface of 3D model.
                  - NODES_Y1: Nodes on NORTH (Y+) surface of 3D model.
                  - NODES_X1: Nodes on EAST (X+) surface of 3D model.
                  - NODES_X0: Nodes on WEST (X-) surface of 3D model.
                  - NODES_Z1: Nodes on TOP (Z+) surface of 3D model.
                  - NODES_Z0: Nodes on BOTTOM (Z-) surface of 3D model.
                These node sets are available for boundary conditions definition.
        * 'ELSET':
                Create boundary element sets. (Default = ON)
                If 'ELSET' is specified, the following element sets are created:
                  - ELEMS_Y0: Elements of SOUTH (Y-) surface of 3D model.
                  - ELEMS_Y1: Elements of NORTH (Y+) surface of 3D model.
                  - ELEMS_X1: Elements of EAST (X+) surface of 3D model.
                  - ELEMS_X0: Elements of WEST (X-) surface of 3D model.
                  - ELEMS_Z1: Elements of TOP (Z+) surface of 3D model.
                  - ELEMS_Z0: Elements of BOTTOM (Z-) surface of 3D model.
        * 'PROPERTY':
                Define an external material mapping law from template file. (Default = None)
                Use in combination with 'matprop' dictionary of material property files and corresponding GV ranges for the material mapping.
    voxelsize : float
        3D model voxelsize.
    eltype : str
        FE element type. The default is eight-node brick element (C3D8 and F3D8). See CalculiX node convention (Par. 6.2.1) at: http://www.dhondt.de/ccx_2.15.pdf
    matpropbits : int
        Bit depth for material mapping.
    refnode
        Reference node coordinates [REF_NODE_x, REF_NODE_y, REF_NODE_z] for kinematic coupling.
        Alternatively use one of the following args [X0, X1, Y0, Y1, Z0, Z1] to generate automatically a REF_NODE at a model boundary.
    """

    # verbose output
    if verbose:
        logging.basicConfig(level=logging.INFO)

    # get 3D dataset shape
    data_shape = voldata.shape

    # voxelsize list if only one value is given
    if type(voxelsize) is not list and type(voxelsize) is not np.ndarray:
        voxelsize = [voxelsize, voxelsize, voxelsize]
    if len(voxelsize) == 1:
        voxelsize = [voxelsize[0], voxelsize[0], voxelsize[0]]

    # variables initialization ##################################################
    # dictionaries (two) of GVs for each user-defined material property
    prop_GV = {}
    prop_GV2 = {}

    # dictionary of the GVs existing in the input model
    existing_GV = []
    GVmin = 0
    GVmax = 2 ** matpropbits - 1

    # dictionary of existing nodes in the model
    nodes = {}
    # dictionary of the model nodal coordinates. For more info see: https://abaqus-docs.mit.edu/2017/English/SIMACAEMODRefMap/simamod-c-node.htm
    nodes_XYZ = {}
    # dictionary of element set nodes
    elset_nodes = {}

    # dictionary of boundary node sets
    nset = {
        'NODES_Y1': [],
        'NODES_Y0': [],
        'NODES_X0': [],
        'NODES_X1': [],
        'NODES_Z1': [],
        'NODES_Z0': []
    }

    refnode_X0 = np.zeros(3)
    refnode_X1 = np.zeros(3)
    refnode_Y0 = np.zeros(3)
    refnode_Y1 = np.zeros(3)
    refnode_Z0 = np.zeros(3)
    refnode_Z1 = np.zeros(3)

    # dictionary of boundary element sets
    elset = {
        'ELEMS_Y1': [],
        'ELEMS_Y0': [],
        'ELEMS_X0': [],
        'ELEMS_X1': [],
        'ELEMS_Z1': [],
        'ELEMS_Z0': []
    }

    # initialize dictionaries of user material properties if given
    if matprop is not None:
        # check matprop dictionary input
        if not "range" in matprop and len(matprop["file"])>1:
            raise IOError('Incorrect material property definition: GV ranges must be given for using multiple property files.')
        elif not "range" in matprop and len(matprop["file"])==1:
            # assign all GVs above zero to the given material property
            matprop["range"] = [[1, GVmax]]

        logging.info('Checking GV range for user material properties')
        n_propfile = 0
        for GVrange in matprop["range"]:
            prop_GVmin, prop_GVmax = int(GVrange[0]), int(GVrange[1])
            # checks on property range Min and Max
            if prop_GVmin < GVmin:
                raise IOError('property GV range below min representable integer')
                # exit(1)

            if prop_GVmax > GVmax:
                raise IOError('property GV range exceeds max representable integer')
                # exit(1)

            prop_GV[n_propfile] = (prop_GVmin, prop_GVmax)
            prop_GV2[n_propfile] = np.arange(prop_GVmin, prop_GVmax + 1, 1)
            n_propfile = n_propfile + 1

        # check prop_GV overlap
        for i in range(n_propfile - 1):
            if not set(prop_GV2[i]).isdisjoint(prop_GV2[i + 1]):
                raise IOError('GV ranges assigned to different material properties cannot overlap')
                # exit(1)

    else:
        if voldata.dtype == 'bool':
            # GV range for binary dataset
            prop_GV[0] = (1, 1)
        else:
            raise IOError('A matprop dictionary must be given for non-binary volume data.')

    # initialize dictionary of elsets and list of materials indexes:
    # for each user-defined material property
    for prop in prop_GV:
        # for each GV assigned to that material property
        for GV in range(prop_GV[prop][0], prop_GV[prop][1] + 1):
            # initialize dictionary of element set nodes corresponding to that GV
            elset_nodes[GV] = []
            # compose list of the existing GVs
            existing_GV.append(GV)

    # setup of element data
    el_i = 0
    col_nodes = data_shape[1] + 1  # n nodes in 1 row
    slice_nodes = (data_shape[2] + 1) * (data_shape[1] + 1)  # n nodes in 1 slice

    logging.info('Composing element sets')
    for slice in range(data_shape[0]):
        for col in range(data_shape[2]):
            for row in range(data_shape[1]):

                # get voxel GV
                GV = voldata[(slice, row, col)]

                if GV in elset_nodes:
                    el_i = el_i + 1

                    # voxel nodes indexes
                    # See eight-node brick element (C3D8 and F3D8) node convention (Par. 6.2.1) at: http://www.dhondt.de/ccx_2.15.pdf
                    el_nodes = [slice_nodes * slice + col_nodes * col + row + 1,
                                slice_nodes * slice + col_nodes * col + row + 2,
                                slice_nodes * slice + col_nodes * (col + 1) + row + 2,
                                slice_nodes * slice + col_nodes * (col + 1) + row + 1,
                                slice_nodes * (slice + 1) + col_nodes * col + row + 1,
                                slice_nodes * (slice + 1) + col_nodes * col + row + 2,
                                slice_nodes * (slice + 1) + col_nodes * (col + 1) + row + 2,
                                slice_nodes * (slice + 1) + col_nodes * (col + 1) + row + 1]

                    # compose dictionary of element sets containing element indexes and element nodes for each set:
                    # for each element set the dictionary contains a list of voxels belonging to that set
                    # the first arg of the list is an element index
                    # the second arg is a list of voxel node indexes
                    elset_nodes[GV].append((el_i, el_nodes))

                    # compose lists of element indexes belonging to boundary sets
                    if slice == 0:
                        elset['ELEMS_Z0'].append(el_i)
                    if slice == data_shape[0] - 1:
                        elset['ELEMS_Z1'].append(el_i)
                    if col == 0:
                        elset['ELEMS_Y0'].append(el_i)
                    if col == data_shape[2] - 1:
                        elset['ELEMS_Y1'].append(el_i)
                    if row == 0:
                        elset['ELEMS_X0'].append(el_i)
                    if row == data_shape[1] - 1:
                        elset['ELEMS_X1'].append(el_i)

                    # store dictionary of existing nodes
                    for i in el_nodes:
                        nodes[i] = 1

    # node coordinates dict; boundary node sets
    logging.info('Detecting node coordinates and boundary nodes')
    node_i = 0
    for slice in range(data_shape[0] + 1):
        for col in range(data_shape[2] + 1):
            for row in range(data_shape[1] + 1):
                node_i = node_i + 1
                if node_i in nodes:
                    # compose dictionary of node coordinates
                    nodes_XYZ[node_i] = (voxelsize[0] * row, voxelsize[1] * col, voxelsize[2] * slice)

                    # compose dictionary of boundary node sets
                    if slice == 0:
                        nset['NODES_Z0'].append(node_i)
                        refnode_Z0 += nodes_XYZ[node_i]
                    if slice == data_shape[0]:
                        nset['NODES_Z1'].append(node_i)
                        refnode_Z1 += nodes_XYZ[node_i]
                    if col == 0:
                        nset['NODES_Y0'].append(node_i)
                        refnode_Y0 += nodes_XYZ[node_i]
                    if col == data_shape[2]:
                        nset['NODES_Y1'].append(node_i)
                        refnode_Y1 += nodes_XYZ[node_i]
                    if row == 0:
                        nset['NODES_X0'].append(node_i)
                        refnode_X0 += nodes_XYZ[node_i]
                    if row == data_shape[1]:
                        nset['NODES_X1'].append(node_i)
                        refnode_X1 += nodes_XYZ[node_i]

    # barycenters of boundary node sets
    refnode_X0 /= len(nset['NODES_X0'])
    refnode_X1 /= len(nset['NODES_X1'])
    refnode_Y0 /= len(nset['NODES_Y0'])
    refnode_Y1 /= len(nset['NODES_Y1'])
    refnode_Z0 /= len(nset['NODES_Z0'])
    refnode_Z1 /= len(nset['NODES_Z1'])

    # # compose dictionary of node coordinates
    # node_i = 0
    # for slice in range(data_shape[0] + 1):
    #     for col in range(data_shape[2] + 1):
    #         for row in range(data_shape[1] + 1):
    #             node_i = node_i + 1
    #             if node_i in nodes:
    #                 nodes_XYZ[node_i] = (voxelsize[0] * row, voxelsize[1] * col, voxelsize[2] * slice)

    # get REF_NODE coordinates
    if refnode:
        if type(refnode) is list:
            if len(refnode) == 1:
                refnode = refnode[0]
            elif len(refnode) != 3:
                logging.warning('Wrong REF_NODE input length. REF_NODE is ignored.')

        if type(refnode) is str:
            if refnode == 'X0':
                refnode = refnode_X0
            elif refnode == 'X1':
                refnode = refnode_X1
            elif refnode == 'Y0':
                refnode = refnode_Y0
            elif refnode == 'Y1':
                refnode = refnode_Y1
            elif refnode == 'Z0':
                refnode = refnode_Z0
            elif refnode == 'Z1':
                refnode = refnode_Z1

    # write ABAQUS *.inp output file ############################################
    logging.info('Start writing INP file')
    # open ABAQUS *.INP output file
    INP = open(fileout, 'w')

    # HEADER
    INP.write('** ---------------------------------------------------------\n')
    INP.write('** Abaqus .INP file written by Voxel_FEM script on {}\n'.format(datetime.now()))
    INP.write('** ---------------------------------------------------------\n')
    # INP.write('*HEADING\n')
    # INP.write('main input {0}\n'.format(filename_in))

    # NODE COORDINATES
    logging.info('Writing model nodes to INP file')
    INP.write('** Node coordinates from input model\n')
    INP.write('*NODE\n')
    node_i = 0

    # write only the existing nodes
    for slice in range(data_shape[0] + 1):
        for col in range(data_shape[2] + 1):
            for row in range(data_shape[1] + 1):
                node_i = node_i + 1
                if node_i in nodes:
                    INP.write('{0:10d}, {n[0]:12.6f}, {n[1]:12.6f}, {n[2]:12.6f}\n'.format(node_i, n=nodes_XYZ[node_i]))

    # ELEMENTS AND ELEMENT SETS
    n_els = 0
    logging.info('Writing model elements to INP file')
    INP.write('** Elements and Element sets from input model\n')
    i = 0

    # for all GVs from template
    for GV in existing_GV:
        i = i + 1

        # if elements with that GV exist
        if len(elset_nodes[GV]) > 0:
            # compose element SET definition string and write to output file
            elset_str = 'SET' + str(GV)
            INP.write('*ELEMENT, TYPE={0}, ELSET={1}\n'.format(eltype, elset_str))

            # write element index followed by list of its nodes
            for el_nodes in elset_nodes[GV]:
                el_i = el_nodes[0]
                INP.write(
                    '{0},{n[0]},{n[1]},{n[2]},{n[3]},{n[4]},{n[5]},{n[6]},{n[7]}\n'.format(el_i, n=el_nodes[1]))
                n_els += 1
                # update dictionary of existing nodes
                for elnd in el_nodes[1]:
                    nodes[elnd] = 1

    # NODE SETS
    if 'NSET' in keywords:
        INP.write('** Additional nset from voxel model. New generated nsets are:\n')
        INP.write('** NODES_Y0, NODES_Y1, NODES_X1, NODES_X0, NODES_Z1, NODES_Z0\n')

        # write node set string
        for nsetName in nset:
            INP.write('*NSET, NSET={0}\n'.format(nsetName))
            CR = 1

            # write node set indexes in lines of 10
            for node_i in nset[nsetName]:
                if CR < 10:
                    INP.write('{},'.format(node_i))
                else:
                    INP.write('{}\n'.format(node_i))
                    CR = 0
                CR = CR + 1
            INP.write('\n')
        logging.info('Additional node sets generated: NODES_Y0, NODES_Y1, NODES_X1, NODES_X0, NODES_Z1, NODES_Z0')

    # ELEMENT SETS
    if 'ELSET' in keywords:
        INP.write('** Additional elset from voxel model. New generated elsets are:\n')
        INP.write('** ELEMS_Y0, ELEMS_Y1, ELEMS_X1, ELEMS_X0, ELEMS_Z1, ELEMS_Z0\n')

        # write element set string
        for elsetName in elset:
            INP.write('*ELSET, ELSET={}\n'.format(elsetName))
            CR = 1

            # write element set indexes in lines of 10
            for el_i in elset[elsetName]:
                if CR < 10:
                    INP.write('{},'.format(el_i))
                else:
                    INP.write('{}\n'.format(el_i))
                    CR = 0
                CR = CR + 1
            INP.write('\n')
        logging.info('Additional element sets generated: ELEMS_Y0, ELEMS_Y1, ELEMS_X1, ELEMS_X0, ELEMS_Z1, ELEMS_Z0')

    # MATERIAL MAPPING
    if 'PROPERTY' in keywords:
        logging.info('User material properties defined. Writing material property section of INP file')
        INP.write('** User material property definition:\n')
        INP.write('** internal variables are: "SetName", "MatName", "GV"\n')

        n_propfile = 0
        for input_property in matprop["file"]:
            # open material property template file
            try:
                PROPfile = open(input_property, 'r')
            except IOError('Material property file not found'):
                exit(1)

            lines = PROPfile.readlines()
            # for each GV corresponding to an existing material index
            for GV in existing_GV:

                # if elements with that GV exist
                if len(elset_nodes[GV]) > 0:

                    # if the GV belongs to the GVrange assigned to the specific material property
                    if prop_GV[n_propfile][0] <= GV <= prop_GV[n_propfile][1]:
                        matset_str = 'SET' + str(GV)

                        for line in lines:
                            line = line.replace('\n', '')

                            # copy comment lines as they are
                            if line.startswith('**'):
                                INP.write('{}\n'.format(line))
                                continue

                            # replace keyword 'SetName' and assign current material to the corresponding (by GV) elset
                            if 'SetName' in line:
                                line = line.replace('SetName', matset_str)

                            # replace keyword 'MatName' and create new material flag for current GV
                            if 'MatName' in line:
                                line = line.replace('MatName', 'MAT' + matset_str)

                            # replace keyword 'GV' with current GV and evaluate expression for each material property
                            if 'GV' in line:
                                prop_strings = line.split(',')
                                prop_n = 0
                                for prop in prop_strings:
                                    if 'GV' in prop:
                                        INP.write('{}'.format(eval(prop)))
                                    else:
                                        INP.write('{}'.format(prop))
                                    if prop_n < len(prop_strings) - 1:
                                        INP.write(',')
                                    prop_n = prop_n + 1

                                INP.write('\n')
                            else:
                                # copy line as is
                                INP.write('{}\n'.format(line))
            PROPfile.close()
            n_propfile = n_propfile + 1

    # BOUNDARY CONDITIONS AND ANALYSIS DEFINITION
    # Open Abaqus analysis template file
    try:
        template = open(templatefile, 'r')
    except IOError('Abaqus template file {} not found.'.format(templatefile)):
        exit(1)
    logging.info('Reading Abaqus template file {}'.format(templatefile))

    # copy line by line info on model solution and boundary conditions from Abaqus analysis template file
    for line in template.readlines():

        if refnode is not None:
            # replace keywords 'refnodeX', 'refnodeY' and 'refnodeZ' with refnode coordinates
            line = line.replace('refnodeX', str(round(refnode[0], 4)))
            line = line.replace('refnodeY', str(round(refnode[1], 4)))
            line = line.replace('refnodeZ', str(round(refnode[2], 4)))

        # write line to output Abaqus file
        INP.write('{}'.format(line))

    template.close()
    INP.close()
    logging.info('Model with {0} nodes and {1} elements written to file {fname}'.format(len(nodes), n_els, fname=fileout))

def matpropdictionary(proplist):
    """Compose dictionary of material properties and property mapping GV ranges.

    Parameters
    ----------
    proplist
        List of material property files followed by the corresponding Gray Value range for material mapping.

    Returns
    -------
    matprop : dict
        Dictionary of material properties for material property mapping:
        matprop = {
            "file": ["prop.inp", "property_temp_bone.inp", ...],
            "range": [[250, 255], [0, 250], ...],
        }
    """

    if proplist is None:
        return None

    if len(proplist) == 0:
        return None

    elif len(proplist) == 1:
        if not os.path.isfile(proplist[0]):
            raise IOError('Property file {0} doesn''t exist.', format(proplist[0]))
        else:
            matprop = {"file": [proplist[0]]}
            return matprop

    else:
        # check correct usage of --mapping flag
        if not len(proplist) % 3 == 0:
            raise IOError('Mapping dictionary error (probably due to incorrect use of --mapping flag): each material property file must be followed by MIN and MAX of the corresponding GV range.')

        else:
            matprop = {
                'file': [],
                'range': []
            }
            for prop in range(0, int(len(proplist)/3)):
                matprop['file'].append(proplist[3*prop])
                matprop['range'].append([int(proplist[3*prop+1]), int(proplist[3*prop+2])])

            return matprop
