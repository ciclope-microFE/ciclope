#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create ABAQUS input file for simulation with the finite element (FE) method from 3D voxel model.

For more information, call this script with the help option::
    stack2abaqus.py -h

"""

"""
2DO:
- jupyter nb - bone microFE example
    - posptrocess CalculiX files
- move TIFF (or other image formats) reading methods to separate library
- documentation pages (sphinx?)
- BC and FE solution parameters definition as command line input (?)

"""

__author__ = 'Gianluca Iori'
__date_created__ = '2020-11-15'
__date__ = '2021-02-14'
__copyright__ = 'Copyright (c) 2021, JC|MSK'
__docformat__ = 'restructuredtext en'
__license__ = "GPL"
__version__ = "0.4.1"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"


import os
import logging
import numpy as np
# import dxchange
import argparse
import textwrap
import tifffile
from datetime import datetime

# parameters #####################################################################


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


def read_tiff_DXchange(filename):
    print('not implemented..')
    # print('Loading file: {}'.format(filename))
    # the dxchange TIFF reader requires list of slice indexes
    # data = dxchange.reader.read_tiff_stack(args.filein, ind=[560, 561, 562, 563, 564])


def main():
    description = textwrap.dedent('''\
            Script to create ABAQUS Finite Element (FE) input file from 3D voxel model.

            The scripts converts the voxels of a 3D input model file to the hexahedra
            of a voxel-FE input file for FE simulations of solid mechanics problems.
            
            The output of this script is an input file (.INP) in ABAQUS syntax
            that can be solved using ABAQUS or CALCULIX.
            
            The script allows the user to define a material mapping strategy
            for the direct conversion of local Grey Values (GVs) of the input 3D model
            to the local material properties of the FE model.
            
            The definition of material mapping laws is done using separate template file(s).
            See "prop.inp" and "property_temp_bone.inp" for examples.

            Boundary conditions (BCs), simulation steps and associated output requests
            are defined in a separate template file. See "tmp.inp" for an example.
            More info at:
            https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-m-Sim-sb.htm#simacae-m-Sim-sb
            ''')

    epilog = textwrap.dedent('''\
            SUPPORTED ABAQUS KEYWORDS:

                For a list of all Abaqus keywords and their description visit:
                https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-c-gen-kwbrowser.htm#simacae-c-gen-kwbrowser__simacae-gen-xsl-U
                
                * 'NODE':
                        Specify nodal coordinates. (Default = Yes)
                * 'ELEMENT':
                        Define elements by giving their nodes. (Default = Yes)
                * 'PROPERTY':
                        Define an external material mapping law from template file. (Default = Yes)
                        Use in combination with '-p' ('--property') and '-r' ('--prange') flags to input
                        material property files and corresponding GV ranges for the material mapping.    
                * 'NSET':
                        Create boundary node sets. (Default = None)
                        If 'NSET' is specified, the following node sets are created:
                          - NODES_S: Nodes on SOUTH surface of 3D model.
                          - NODES_N: Nodes on NORTH surface of 3D model.
                          - NODES_E: Nodes on EAST surface of 3D model.
                          - NODES_W: Nodes on WEST surface of 3D model.
                          - NODES_T: Nodes on TOP surface of 3D model.
                          - NODES_B: Nodes on BOTTOM surface of 3D model.
                        These node sets are available for boundary conditions definition.
                * 'ELSET':
                        Create boundary element sets. (Default = None)
                        If 'ELSET' is specified, the following element sets are created:
                          - ELEMS_S: Elements of SOUTH surface of 3D model.
                          - ELEMS_N: Elements of NORTH surface of 3D model.
                          - ELEMS_E: Elements of EAST surface of 3D model.
                          - ELEMS_W: Elements of WEST surface of 3D model.
                          - ELEMS_T: Elements of TOP surface of 3D model.
                          - ELEMS_B: Elements of BOTTOM surface of 3D model.
                
            EXAMPLES:
            
            * Convert CT scan of embedded bone tissue mapping the model GVs to local material properties.
              GVs between 251 and 255 are reserved for the embedding matrix, and a separate property file is given
              for the mechanical properties of this material. The type of FE analysis, all BCs and requested output are
              defined in the template file "tmp.inp"::    
                
                stack2abaqus.py scan-sample_15-3-crop_c0001.tif pippo.inp
                -k NODE ELEMENT NSET PROPERTY
                -p ./material_properties/bone.inp ./material_properties/PMMA.inp
                -pr 1:250 251:255
                -t input_templates/tmp.inp
            ''')

    parser = argparse.ArgumentParser(description = description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filein', type=str, help='Input filename (voxel data).')
    parser.add_argument('fileout', type=str, default=None, help='Output filename (Abaqus .INP).')
    parser.add_argument('-vs', '--voxelsize', type=float, nargs='+', default=[1.0, 1.0, 1.0], help='Data voxel size [X, Y, Z].')
    parser.add_argument('-k', '--keywords', type=str, nargs='*', default=['NODE', 'ELEMENT', 'PROPERTY'], help='Abaqus keywords.')
    parser.add_argument('--bits', type=int, default=8, help='Bit depth for material mapping.')
    parser.add_argument('--eltype', type=str, default='C3D8', help='Element type.')
    parser.add_argument('-p', '--property', type=str, nargs='*', default='prop.inp', help='Template file for material property mapping.')
    parser.add_argument('-pr', '--prange', type=str, nargs='*', default='1:255', help='GV range for user material property.')
    parser.add_argument('-t', '--template', type=str, default='input_templates/tmp.inp', help='Template file (Abaqus syntax) defining analysis steps, boundary conditions and output requests.')
    parser.add_argument('-v', '--verbose', type=bool, default=False, help='Verbose output')

    args = parser.parse_args()

    # check inputs
    filename_in, ext = os.path.split(args.filein)
    if args.verbose is True:
        logging.basicConfig(level=logging.INFO)

    # check if property and prange have the same number of inputs
    if args.property is not None and args.prange is not None:
        if len(args.property) != len(args.prange):
            raise IOError('The length of property and prange must be the same')
            # exit(1)

    # load 3D voxel data
    data = read_tiff_stack(args.filein)
    data_shape = data.shape
    logging.info('Data loaded with size {size[0]} x {size[1]} x {size[2]}'.format(size=data_shape))

    # variables initialization
    prop_GV = {} # dictionaries (two) of GVs for each user-defined material property
    prop_GV2 = {}
    existing_GV = [] # dictionary of the GVs existing in the input model
    GVmin = 0
    GVmax = 2**args.bits - 1
    nodes = {} # dictionary of existing nodes of the output FE model
    nodes_XYZ = {}  # dictionary of the model nodal coordinates. For more info see: https://abaqus-docs.mit.edu/2017/English/SIMACAEMODRefMap/simamod-c-node.htm
    elset_nodes = {} # dictionary of element set nodes

    # dictionary of boundary node sets
    nset = {
      'NODES_N': [],
      'NODES_S': [],
      'NODES_W': [],
      'NODES_E': [],
      'NODES_T': [],
      'NODES_B': []
    }

    # dictionary of boundary element sets
    elset = {
      'ELEMS_N': [],
      'ELEMS_S': [],
      'ELEMS_W': [],
      'ELEMS_E': [],
      'ELEMS_T': [],
      'ELEMS_B': []
    }

    # initialize dictionaries of user material properties if given
    if args.property is not None:
        n_propfile = 0
        for prange in args.prange:
            GVrange = prange.split(':')
            prop_GVmin, prop_GVmax = int(GVrange[0]), int(GVrange[1])
            # checks on property range Min and Max
            if prop_GVmin < GVmin:
                raise IOError('property GV range below representable grey values')
                # exit(1)

            if prop_GVmax > GVmax:
                raise IOError('property GV range exceeds max. representable grey value')
                # exit(1)

            prop_GV[n_propfile] = (prop_GVmin, prop_GVmax)
            prop_GV2[n_propfile] = np.arange(prop_GVmin, prop_GVmax + 1, 1)
            n_propfile = n_propfile + 1

    # check prop_GV overlap
    for i in range(n_propfile-1):
        if not set(prop_GV2[i]).isdisjoint(prop_GV2[i+1]):
            raise IOError('GV ranges assigned to different material properties cannot overlap')
            # exit(1)

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
    col_nodes = data_shape[1] + 1 # n nodes in 1 row
    slice_nodes = (data_shape[2] + 1) * (data_shape[1] + 1) # n nodes in 1 slice

    for slice in range(data_shape[0]):
        for col in range(data_shape[2]):
            for row in range(data_shape[1]):

                # get voxel GV
                GV = data[(slice, col, row)]
                
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
                        elset['ELEMS_B'].append(el_i)
                    if slice == data_shape[0] - 1:
                        elset['ELEMS_T'].append(el_i)
                    if col == 0:
                        elset['ELEMS_S'].append(el_i)
                    if col == data_shape[2] - 1:
                        elset['ELEMS_N'].append(el_i)
                    if row == 0:
                        elset['ELEMS_W'].append(el_i)
                    if row == data_shape[1] - 1:
                        elset['ELEMS_E'].append(el_i)

                    # store dictionary of existing nodes
                    for i in el_nodes:
                        nodes[i] = 1

    # compose lists of node indexes belonging to boundary sets
    node_i = 0
    for slice in range(data_shape[0] + 1):
        for col in range(data_shape[2] + 1):
            for row in range(data_shape[1] + 1):
                node_i = node_i + 1
                if node_i in nodes:
                    if slice == 0:
                        nset['NODES_B'].append(node_i)
                    if slice == data_shape[0]:
                        nset['NODES_T'].append(node_i)
                    if col == 0:
                        nset['NODES_S'].append(node_i)
                    if col == data_shape[2]:
                        nset['NODES_N'].append(node_i)
                    if row == 0:
                        nset['NODES_W'].append(node_i)
                    if row == data_shape[1]:
                        nset['NODES_E'].append(node_i)

    # compose dictionary of node coordinates
    node_i = 0
    for slice in range(data_shape[0] + 1):
        for col in range(data_shape[2] + 1):
            for row in range(data_shape[1] + 1):
                node_i = node_i + 1
                if node_i in nodes:
                    nodes_XYZ[node_i] = (args.voxelsize[0] * row, args.voxelsize[1] * col, args.voxelsize[2] * slice)

    # open ABAQUS *.INP output file
    INP = open(args.fileout, 'w')

    # write ABAQUS *.inp output file
    # HEADER:
    INP.write('** ---------------------------------------------------------\n')
    INP.write('** Abaqus .INP file written on {}\n'.format(datetime.now()))
    INP.write('** ---------------------------------------------------------\n')
    INP.write('*HEADING\n')
    INP.write('main input {0}\n'.format(filename_in))

    # NODAL COORDINATES:
    if 'NODE' in args.keywords:
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

    # ELEMENTS AND ELEMENT SETS:
    if 'ELEMENT' in args.keywords:
        INP.write('** Elements and Element sets from input model\n')
        i = 0

        # for all GVs from template
        for GV in existing_GV:
            i = i + 1

            # if elements with that GV exist
            if len(elset_nodes[GV]) > 0:
                # compose element SET definition string and write to output file
                elset_str = 'SET' + str(GV)
                INP.write('*ELEMENT, TYPE={0}, ELSET={1}\n'.format(args.eltype, elset_str))

                # write element index followed by list of its nodes
                for el_nodes in elset_nodes[GV]:
                    el_i = el_nodes[0]
                    INP.write('{0},{n[0]},{n[1]},{n[2]},{n[3]},{n[4]},{n[5]},{n[6]},{n[7]}\n'.format(el_i, n=el_nodes[1]))
                    # update dictionary of existing nodes
                    for elnd in el_nodes[1]:
                        nodes[elnd] = 1

    # NODE SETS:
    if 'NSET' in args.keywords:
        INP.write('** Additional nset from voxel model. New generated nsets are:\n')
        INP.write('** NODES_S, NODES_N, NODES_E, NODES_W, NODES_T, NODES_B\n')

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

    # ELEMENT SETS:
    if 'ELSET' in args.keywords:
        INP.write('** Additional elset from voxel model. New generated elsets are:\n')
        INP.write('** ELEMS_S, ELEMS_N, ELEMS_E, ELEMS_W, ELEMS_T, ELEMS_B\n')

        # write element set string
        for elsetName in elset:
            INP.write('*ELSET, ELSET={}\n'.format(elsetName))
            CR = 1

            # write element set indexes in lines of 10
            for el_i in elset[elsetName]:
                if CR < 10:
                    INP.write('{}'.format(el_i))
                else:
                    INP.write('{}\n'.format(el_i))
                    CR = 0
                CR = CR + 1
            INP.write('\n')

    # MATERIAL MAPPING:
    if 'PROPERTY' in args.keywords:
        logging.info('User material properties defined.')
        INP.write('** User material property definition:\n')
        INP.write('** internal variables are: "SetName", "MatName", "GV"\n')

        n_propfile = 0
        for property in args.property:
            # open material property template file
            try:
                PROPfile = open(property, 'r')
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
                                INP.write('{}\n'.format(line))
            PROPfile.close()
            n_propfile = n_propfile + 1

    # copy line by line info on model solution and boundary conditions from Abaqus template file
    try:
        BCfile = open(args.template, 'r')
    except IOError('Abaqus template file {} not found.'.format(args.template)):
        exit(1)
    logging.info('Reading Abaqus template file {}'.format(args.template))

    for line in BCfile.readlines():
        # copy line to output Abaqus file
        INP.write('{}'.format(line))

    BCfile.close()
    INP.close()
    logging.info('Data written to file {}'.format(args.fileout))
    return

if __name__ == '__main__':
    main()
