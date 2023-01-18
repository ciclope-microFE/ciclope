#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
From 3D Computed Tomography (CT) images to Finite Element (FE) models.
For more information, call this script with the help option:
    ciclope.py -h

"""

__author__ = ['Gianluca Iori', 'Martino Pani']
__date_created__ = '2021-08-06'
__date__ = '2023-01-18'
__copyright__ = 'Copyright (c) 2023, ORMIR'
__docformat__ = 'restructuredtext en'
__license__ = "MIT"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"

import os
import argparse
import logging
import textwrap
import numpy as np
import meshio
from skimage.filters import gaussian
from ciclope.core import voxelFE, tetraFE
from ciclope.utils import preprocess, recon_utils
from sys import version_info
from .__about__ import __version__

#################################################################################

def main():
    description = textwrap.dedent('''\
                Ciclope generates CalculiX Finite Element (FE) input files from 3D voxel data.

                The output of this script is an input file (.INP) in ABAQUS syntax
                that can be solved using ABAQUS or CALCULIX.

                The input 3D image is segmented and processed to maintain only the largest connected structure.
                3D voxel or tetrahedra FE input files can be generated.
                Meshes for 3D rendering of the outer model shell and of the whole selected volume can be generated.

                The user can define a material mapping strategy (implemented only for voxel FE)
                for the direct conversion of local Grey Values (GVs) of the 3D image
                to local material properties of the FE model.

                The definition of material mapping laws is done using separate template file(s).
                See the folder "material_properties" for examples.

                Boundary conditions (BCs), analysis steps and associated output requests
                are defined in a separate template file. See the folder "input_templates" for examples.
                More info on how to define analysis steps in ABAQUS at:
                https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-m-Sim-sb.htm#simacae-m-Sim-sb
                ''')
    epilog = textwrap.dedent('''\
                EXAMPLES:

                * Convert SRCT scan of single strut to voxel FE model.
                  GVs between 1 and 249 are used for mapping local bone tissue material properties.
                  GVs between 250 and 255 are reserved for a steel cap, and a separate property file is given
                  for the mechanical properties of steel. A static tensile analysis and output requests are
                  defined in the ABAQUS template file "tmp_example02_tens_static.inp":

                    ciclope.py /home/gianthk/Data/TOMCAT/Kaya/D_single_h1h2_scale05/D_single_h1h2_scale050017.tif pippo.inp
                    -vs 3.25 3.25 3.25
                    -r 5
                    --voxelfe
                    --template /home/gianthk/PycharmProjects/CT2FE/input_templates/tmp_example01_tens_static.inp
                    --mapping
                    /home/gianthk/PycharmProjects/CT2FE/material_properties/steel.inp 250 255
                    /home/gianthk/PycharmProjects/CT2FE/material_properties/bone.inp 1 249

                * Same as the example above, but without material mapping.
                  The image is segmented and all voxels belonging to the largest connected strut are assigned
                  to the same elements set with uniform material properties. The analysis template file 
                   is assumed to contain the material property definition (see file "tmp_example02_tens_static.inp"):

                    ciclope.py /home/gianthk/Data/TOMCAT/Kaya/D_single_h1h2_scale05/D_single_h1h2_scale050017.tif pippo.inp
                    -vs 3.25 3.25 3.25
                    -r 5
                    --voxelfe
                    --template /home/gianthk/PycharmProjects/CT2FE/input_templates/tmp_example02_tens_static.inp
                ''')
    onemoreexample = textwrap.dedent('''
                * Convert CT scan of embedded bone tissue mapping the model GVs to local material properties.
                  GVs between 251 and 255 are reserved for the embedding matrix, and a separate property file is given
                  for the mechanical properties of this material. The type of FE analysis, all BCs and requested output are
                  defined in the template file "tmp.inp":

                    ciclope.py scan-sample_15-3-crop_c0001.tif pippo.inp
                    -r 4
                    --vol_mesh
                    --shell_mesh
                    --voxelfe
                    --template /home/gianthk/PycharmProjects/CT2FE/input_templates/tmp_example02_tens_static.inp
                    ''')

    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filein', type=str, help='<Required> Input filename (voxel data).')
    parser.add_argument('fileout', type=str, default=None, help='<Required> Output filename (Abaqus .INP).')
    parser.add_argument('-vs', '--voxelsize', type=float, default=[1., 1., 1.], nargs='+', help='Voxel size.')
    parser.add_argument('-r', '--resampling', type=float, default=1., help='Resampling factor.')
    parser.add_argument('-t', '--threshold', type=int, default=None, help='Threshold value.')
    parser.add_argument('-s', '--smooth', type=float, nargs='?', const=1., default=0.,
                        help='Smooth image with gaussian filter of given Sigma before thresholding.')
    parser.add_argument('--caps', type=int, default=None,
                        help='Add caps of given thickness to the bottom and top of the model for mesh creation.')
    parser.add_argument('--caps_val', type=int, default=0, help='Caps grey value.')
    parser.add_argument('--shell_mesh', dest='shell_mesh', action='store_true',
                        help='Write VTK mesh of outer shell generated with PyMCubes.')
    parser.add_argument('--vol_mesh', dest='vol_mesh', action='store_true',
                        help='Write VTK volume mesh of tetrahedra with pygalmesh.')
    parser.add_argument('--max_facet_distance', type=float, default=None, help='CGAL mesh parameter.')
    parser.add_argument('--max_cell_circumradius', type=float, default=None, help='CGAL mesh parameter.')
    parser.add_argument('--voxelfe', dest='voxelfe', action='store_true', help='Write voxel FE model (.INP) file.')
    parser.add_argument('--template', type=str, default=None,
                        help='<Required by --voxelfe> Abaqus analysis template file (.INP).')
    parser.add_argument('-m', '--mapping', default=None, nargs='+',
                        help='Template file for material property mapping. If more than one property is given, each property filename must followed by the corresponding GV range.')
    parser.add_argument('--tetrafe', dest='tetrafe', action='store_true',
                        help='Write linear tetrahedra FE model (.INP) file.')
    parser.add_argument('--refnode', default=None, nargs='+',
                        help='Reference node input. Used for kinematic coupling of Boundary Conditions in the analysis template file.'
                             'The REF_NODE coordinates [x,y,z] can be given. Alternatively use one of the following args [X0, X1, Y0, Y1, Z0, Z1] to generate a REF_NODE at a model boundary.')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose output.')
    parser.add_argument('--version', action='version', version=_get_version_text(), help='Display version information.')
    parser.set_defaults(shell_mesh=False, vol_mesh=False, voxelfe=False, tetrafe=False, verbose=False)

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    # filename base
    [fileout_base, type] = os.path.splitext(args.fileout)

    # isotropic voxel size (spacing) [mm]
    sp = np.array(args.voxelsize)

    # Read tiff stack #######################################################
    I = recon_utils.read_tiff_stack(args.filein)

    # Gaussian smooth #######################################################
    if args.smooth != 0:
        logging.info('Smoothing image with gaussian kernel')
        I = gaussian(I, sigma=args.smooth, preserve_range=True)

    # Resize the dataset ####################################################
    # use scikit
    if args.resampling != 1:
        logging.info("Resampling input dataset with factor 1/{}".format(args.resampling))
        I, sp = preprocess.resample(I, sp, args.resampling)
        logging.info("New voxelsize: {}".format(sp))

    if args.verbose:
        # plot the midplanes through the stack
        recon_utils.plot_midplanes(I)
        # write midplanes as .PNG
        recon_utils.writemidplanes(I, fileout_base + ".png")

    # Add image caps before thresholding ####################################
    if args.caps is not None:
        # add cap so that the mesh is closed
        I = recon_utils.add_cap(I, cap_thickness=args.caps, cap_val=args.caps_val)

    # Binarise the dataset ##################################################
    BW, T = preprocess.segment(I, args.threshold)
    if args.threshold is None:
        logging.info("Threshold value not given. Using Otsu method..")

    logging.info("Threshold: {}".format(T))

    # Keep largest isolated cluster of voxels #################################
    logging.info("Removing unconnected clusters of voxels..")
    L = preprocess.remove_unconnected(BW)

    # Visualization with Napari #############################################
    # import napari
    # viewer = napari.view_image(L)

    # Generate outer shell mesh #############################################
    if args.shell_mesh:
        logging.info("Writing triangle mesh of outer shell")
        # shell mesh using pymcubes (high resolution mesh; caps excluded)
        vertices, triangles, shellmesh = tetraFE.shell_mesh(L, method='pymcubes')

        # write VTK mesh with meshio
        meshio.write_points_cells(fileout_base + "_shell.vtk", vertices.tolist(), [("triangle", triangles.tolist())])

        # shell mesh using pygalmesh (CGAL) (low resolution mesh; caps included)
        # check CGAL parameters
        # max_facet_distance, max_cell_circumradius = check_cgal_params(args.max_facet_distance, args.max_cell_circumradius, sp)

        # vertices, triangles, shellmesh = shell_mesh(L, 'pygalmesh', sp, max_facet_distance, max_cell_circumradius)
        # shellmesh.write(fileout_base + "_shell.vtk")

    # Generate volume mesh ##################################################
    if args.vol_mesh:
        # check CGAL parameters
        max_facet_distance, max_cell_circumradius = tetraFE.check_cgal_params(args.max_facet_distance,
                                                                              args.max_cell_circumradius, sp)

        # generate mesh with CGAL
        volmesh = tetraFE.cgal_mesh(L, sp, 'tetra', max_facet_distance, max_cell_circumradius)

        # write the mesh to file
        logging.info("Writing tetrahedra volume mesh")
        volmesh.write(fileout_base + "_tetramesh.vtk")

    # Generate voxel FE model ##################################################
    if args.voxelfe:
        # check voxelfe args
        if not args.template:
            raise IOError('Analysis template file (--template) required with --voxelfe.')

        # build dictionary of material properties
        matprop = voxelFE.matpropdictionary(args.mapping)

        # generate unstructured grid mesh from 3D data
        mesh, refnodes = voxelFE.vol2ugrid(L, sp, refnodes=True)

        if args.refnode is not None:
            if args.refnode[0] in refnodes:
                refnode = refnodes[args.refnode[0]]
            elif len(args.refnode) > 1:
                refnode = args.refnode
            else:
                logging.warning('Refnode not recognized. Option deactivated.')
                refnode = None
        else:
            refnode = args.refnode

        if args.mapping is None:
            # voxelFE of binary volume data; the material property definition is assumed to be in the analysis template file
            voxelFE.mesh2voxelfe(mesh, args.template, fileout_base + "_voxelFE.inp", keywords=['NSET', 'ELSET'], refnode=refnode, verbose=args.verbose)

        else:
            # mask the data (keep only the largest connected volume)
            masked = I
            masked[~L] = 0
            # voxelFE of greyscale volume data; material mapping on
            voxelFE.mesh2voxelfe(mesh, args.template, fileout_base + "_voxelFE.inp", matprop, keywords=['NSET', 'ELSET', 'PROPERTY'], refnode=refnode, verbose=args.verbose)

    # Generate tetrahedra FE model #############################################
    if args.tetrafe:
        # check tetrafe args
        if not args.template:
            raise IOError('Analysis template file (--template) required with --tetrafe.')

        # # build dictionary of material properties
        # matprop = matpropdictionary(args.mapping)

        if args.mapping is None:
            if 'volmesh' not in locals():
                # check CGAL parameters
                max_facet_distance, max_cell_circumradius = tetraFE.check_cgal_params(args.max_facet_distance, args.max_cell_circumradius, sp)

                # generate volume mesh
                volmesh = tetraFE.cgal_mesh(L, sp, 'tetra', max_facet_distance, max_cell_circumradius)

            elif volmesh is None:
                # check CGAL parameters
                max_facet_distance, max_cell_circumradius = tetraFE.check_cgal_params(args.max_facet_distance, args.max_cell_circumradius, sp)

                # generate volume mesh
                volmesh = tetraFE.cgal_mesh(L, sp, 'tetra', max_facet_distance, max_cell_circumradius)

            # tetraFE of the already meshed volume; the material property definition is assumed to be in the analysis template file
            tetraFE.mesh2tetrafe(volmesh, args.template, fileout_base + "_tetraFE.inp", keywords=['NSET', 'ELSET'],
                                 verbose=args.verbose)

        else:
            logging.exception("--tetrafe with mapping option not implemented yet.")

    return

def _get_version_text():
    python_version = f"{version_info.major}.{version_info.minor}.{version_info.micro}"
    return "\n".join(
        [
            f"ciclope {__version__} [Python {python_version}]",
            "Copyright (c) 2023, Gianluca Iori et al.",
        ]
    )

if __name__ == '__main__':
    main()



