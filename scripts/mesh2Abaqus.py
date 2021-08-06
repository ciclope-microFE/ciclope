#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create ABAQUS input file for simulation with the finite element (FE) method from unstructured grid VTK mesh.

For more information, call this script with the help option::
    mesh2abaqus.py -h

"""

"""
2DO:
- jupyter notebook documenting one full example
- run complete test with Calculix
- documentation pages (sphinx?)
- BC and FE solution parameters definition as command line input (?)

"""

__author__ = 'Gianluca Iori'
__date_created__ = '2021-01-13'
__date__ = '2021-01-14'
__copyright__ = 'Copyright (c) 2021, JC|MSK'
__docformat__ = 'restructuredtext en'
__license__ = "GPL"
__version__ = "0.4.1"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"


import logging
import numpy as np
import argparse
import textwrap
import meshio

# parameters #####################################################################


def main():
    description = textwrap.dedent('''\
            Script to create ABAQUS Finite Element (FE) input file from 3D unstructured grid VTK mesh.

            The scripts converts a .VTK unstructured grid mesh file to an Abaqus .INP file.            
            The output of this script is an input file (.INP) in ABAQUS syntax
            that can be solved using ABAQUS or CalculiX.
            
            Boundary conditions (BCs), simulation steps and associated output requests
            are defined in a separate template file. See "tmp.inp" for an example.
            More info at:
            https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-m-Sim-sb.htm#simacae-m-Sim-sb
            ''')

    epilog = textwrap.dedent('''\
            SUPPORTED ABAQUS KEYWORDS:

                For a list of all Abaqus keywords and their description visit:
                https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-c-gen-kwbrowser.htm#simacae-c-gen-kwbrowser__simacae-gen-xsl-U
                
                * 'PROPERTY' (Not supported):
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
            
            * Convert VTK unstructured grid mesh to CalculiX FE for tensile test simulations.
              The analysis type and all BCs and requested outputs are defined in the template file "tmp_example02_tens_Nlgeom.inp":
                
                mesh2Abaqus.py D_single.vtk CalculiX/D_single.inp
                -k NSET ELSET
                -t input_templates/tmp_example02_tens_Nlgeom.inp

            ''')

    parser = argparse.ArgumentParser(description = description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filein', type=str, help='Input filename (voxel data).')
    parser.add_argument('fileout', type=str, default=None, help='Output filename (Abaqus .INP).')
    parser.add_argument('-k', '--keywords', type=str, nargs='*', default=None, help='Abaqus keywords.')
    # parser.add_argument('--eltype', type=str, default='C3D8', help='Element type.')
    parser.add_argument('-t', '--template', type=str, default='input_templates/tmp.inp', help='Template file (Abaqus syntax) defining analysis steps, boundary conditions and output requests.')
    parser.add_argument('--float_fmt', type=str, default='.6e', help='Float formatting option.')
    parser.add_argument('-v', '--verbose', type=bool, default=False, help='Verbose output')

    args = parser.parse_args()

    # verbose output
    if args.verbose is True:
        logging.basicConfig(level=logging.INFO)

    # load 3D mesh using meshio module
    mesh = meshio.read(args.filein)
    logging.info('Loaded mesh with fields:\n')
    logging.info(vars(mesh))

    if 'NSET' in args.keywords:
        # find model boundaries
        model_coors_max = np.amax(mesh.points, 0)  # [microns?]
        model_coors_min = np.amin(mesh.points, 0)

        # add dictionary of boundary point sets
        mesh.point_sets = {
            'NODES_N': np.where(mesh.points[:, 1] == model_coors_max[1])[0],
            'NODES_S': np.where(mesh.points[:, 1] == model_coors_min[1])[0],
            'NODES_W': np.where(mesh.points[:, 0] == model_coors_max[0])[0],
            'NODES_E': np.where(mesh.points[:, 0] == model_coors_min[0])[0],
            'NODES_T': np.where(mesh.points[:, 2] == model_coors_max[2])[0],
            'NODES_B': np.where(mesh.points[:, 2] == model_coors_min[2])[0]
        }

    if 'ELSET' in args.keywords:
        # dictionary of boundary element sets
        elset = {
          'ELEMS_N': [],
          'ELEMS_S': [],
          'ELEMS_W': [],
          'ELEMS_E': [],
          'ELEMS_T': [],
          'ELEMS_B': []
        }

        # find boundary cells
        ELEMS_E = np.array([]).astype('int')
        ELEMS_W = np.array([]).astype('int')
        ELEMS_S = np.array([]).astype('int')
        ELEMS_N = np.array([]).astype('int')
        ELEMS_T = np.array([]).astype('int')
        ELEMS_B = np.array([]).astype('int')

        for node_e in mesh.point_sets['NODES_E']:
            ELEMS_E = np.append(ELEMS_E, np.where(np.any(mesh.cells[0][1] == node_e, axis=1)))

        for node_w in mesh.point_sets['NODES_W']:
            ELEMS_W = np.append(ELEMS_W, np.where(np.any(mesh.cells[0][1] == node_w, axis=1)))

        for node_s in mesh.point_sets['NODES_S']:
            ELEMS_S = np.append(ELEMS_S, np.where(np.any(mesh.cells[0][1] == node_s, axis=1)))

        for node_n in mesh.point_sets['NODES_N']:
            ELEMS_N = np.append(ELEMS_N, np.where(np.any(mesh.cells[0][1] == node_n, axis=1)))

        for node_t in mesh.point_sets['NODES_T']:
            ELEMS_T = np.append(ELEMS_T, np.where(np.any(mesh.cells[0][1] == node_t, axis=1)))

        for node_b in mesh.point_sets['NODES_B']:
            ELEMS_B = np.append(ELEMS_B, np.where(np.any(mesh.cells[0][1] == node_b, axis=1)))

        mesh.cell_sets = {
            'ELEMS_N': [np.unique(ELEMS_N)],
            'ELEMS_S': [np.unique(ELEMS_S)],
            'ELEMS_W': [np.unique(ELEMS_W)],
            'ELEMS_E': [np.unique(ELEMS_E)],
            'ELEMS_T': [np.unique(ELEMS_T)],
            'ELEMS_B': [np.unique(ELEMS_B)],
            'SET1': [np.arange(0, len(mesh.cells[0][1]))]
        }

    else:
        mesh.cell_sets = {'SET1': [np.arange(0, len(mesh.cells[0][1]))]}

    # write Abaqus mesh using meshio
    meshio.abaqus.write(args.fileout, mesh, args.float_fmt)

    # Add analysis definition to the end of the Abaqus INP from template
    # Open ABAQUS *.INP output file in Append mode
    INP = open(args.fileout, 'a')

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
