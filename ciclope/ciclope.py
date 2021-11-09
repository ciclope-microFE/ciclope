#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
From 3D Computed Tomography (CT) images to Finite Element (FE) models.
For more information, call this script with the help option:
    ciclope.py -h

"""

__author__ = ['Gianluca Iori', 'Martino Pani']
__date_created__ = '2021-08-06'
__date__ = '2021-11-03'
__copyright__ = 'Copyright (c) 2021, JC|MSK'
__docformat__ = 'restructuredtext en'
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = 'Gianluca Iori'
__email__ = "gianthk.iori@gmail.com"

import os
import argparse
import logging
import textwrap
import numpy as np
import recon_utils as ru
import meshio
import mcubes
from scipy import ndimage, misc
from skimage.filters import threshold_otsu, gaussian
# from skimage import measure, morphology
from pybonemorph import remove_unconnected
import matplotlib.pyplot as plt
from datetime import datetime

#################################################################################

def shell_mesh(bwimage, method='pymcubes', voxelsize=[1., 1., 1.], max_facet_distance=0.0, max_cell_circumradius=0.0):
    """Generate outer shell mesh of triangles from binary volume data.
    The mesh is generated using the PyMCubes module and the smooth function contained in it:
    https://github.com/pmneila/PyMCubes

    Alternatively, the marching cube algorithm from the scikit-image python module can be used:
    https://scikit-image.org/docs/dev/api/skimage.measure.html?highlight=marching#skimage.measure.marching_cubes


    Parameters
    ----------
    bwimage
        Binary image.
    method : str
        'pymcubes': PyMCubes module.
        'marching_cubes': scikit-image's marching cube algorithm.
        'pygalmesh': pygalmesh module (CGAL).
    voxelsize : float
        Image voxelsize.
    max_facet_distance : float
        CGAL parameter.
    max_cell_circumradius : float
        CGAL parameter.

    Returns
    -------
    vertices
        Mesh vertices.
    triangles
        Mesh triangles.
    shellmesh : meshio
        meshio mesh
    """

    if method == 'pymcubes':
        # (the 0-levelset of the output of mcubes.smooth is the smoothed version of the 0.5-levelset of the binary array.
        smoothed_L = mcubes.smooth(bwimage)
        # Extract the 0-levelset
        vertices, triangles = mcubes.marching_cubes(np.transpose(smoothed_L, [2, 1, 0]), 0)
        shellmesh = None

    elif method == 'marching_cubes':
        from skimage.measure import marching_cubes
        vertices, triangles, tmp, tmp2 = marching_cubes(np.transpose(bwimage, [2, 1, 0]), level=None, step_size=1)
        shellmesh = None

    elif method == 'pygalmesh':
        shellmesh = cgal_mesh(bwimage, voxelsize, 'triangle', max_facet_distance, max_cell_circumradius)
        vertices = None
        triangles = None

    else:
        raise IOError('{0} method unknown.', format(method))

    # get SOUTH cap vertices
    # import pygmsh
    # vertices_capS = vertices[vertices[:, 0] == vertices[:, 0].min()]
    # with pygmsh.geo.Geometry() as geom:
    #     geom.add_polygon(vertices_capS[:,1:3], mesh_size=1)
    #     mesh_capS = geom.generate_mesh()
    #
    # import pygalmesh
    # mesh_capS = pygalmesh.generate_2d(vertices_capS[:,1:3],max_edge_size=0.2,num_lloyd_steps=10)

    return vertices, triangles, shellmesh

def check_cgal_params(max_facet_distance, max_cell_circumradius, voxelsize):
    """Check CGAL mesher parameters.
    # https://github.com/nschloe/pygalmesh#volume-meshes-from-surface-meshes

    Parameters
    ----------
    max_facet_distance
    max_cell_circumradius

    voxelsize : float
        Image voxel size.

    Returns
    -------
    max_facet_distance : float
    max_cell_circumradius : float
    """

    if len(voxelsize) > 1:
        voxelsize = voxelsize[0]

    if max_facet_distance is None:
        max_facet_distance = 1 * voxelsize

    if max_cell_circumradius is None:
        max_cell_circumradius = 5 * voxelsize

    return max_facet_distance, max_cell_circumradius

def cgal_mesh(bwimage, voxelsize, meshtype='both', max_facet_distance=0.0, max_cell_circumradius=0.0):
    """Generate mesh of from binary volume data using CGAL.
    The mesh is generated using the PyGalmesh module. For more info visit: https://github.com/nschloe/pygalmesh#volume-meshes-from-surface-meshes
    The pygalmesh.generate_from_array method returns a mesh containing both a cells set of tetrahedra (volume mesh)
    and a cells set of triangles (shell mesh). The parameter 'meshtype' is used to control which type of mesh is returned.

    Parameters
    ----------
    bwimage
        Binary image.
    voxelsize : float
        Image voxelsize.
    meshtype : str
        'triangle': Outer mesh (shell) of triangles.
        'tetra': Volume mesh of tetrahedra.
        'both': Both shell and volume cells sets.
    max_facet_distance : float
        CGAL parameter.
    max_cell_circumradius : float
        CGAL parameter.

    Returns
    -------
    mesh

    """

    import pygalmesh
    # generate mesh
    mesh = pygalmesh.generate_from_array(np.transpose(bwimage, [2, 1, 0]).astype('uint8'), tuple(voxelsize), max_facet_distance=max_facet_distance, max_cell_circumradius=max_cell_circumradius)

    # check mesh fields and get index for triangle and tetra cells
    triangle_id = None
    tetra_id = None

    for id in range(0, len(mesh.cells)):
        if mesh.cells[id].type == 'triangle':
            triangle_id = id
        elif mesh.cells[id].type == 'tetra':
            tetra_id = id

    # pygalmesh returns two cells entities ('triangle' and 'tetra'). Select one or both for return
    if len(mesh.cells) > 1:
        if meshtype == 'triangle':
            if triangle_id is None:
                raise Exception('pygalmesh did not return a triangle mesh.')
            else:
                logging.info("Removing tetra cells.")
                mesh.cells = [mesh.cells[triangle_id]]
                # this is prone to errors..
                mesh.cell_data['medit:ref'] = mesh.cell_data['medit:ref'][triangle_id]
        elif meshtype == 'tetra':
            if tetra_id is None:
                raise Exception('pygalmesh did not return a tetra mesh.')
            else:
                logging.info("Removing triangle cells.")
                mesh.cells = [mesh.cells[tetra_id]]
                # this is prone to errors..
                mesh.cell_data['medit:ref'] = mesh.cell_data['medit:ref'][tetra_id]
        elif meshtype == 'both':
            logging.info("Both triangle and tetra cells kept.")
        else:
            raise IOError('{0} method type.', format(meshtype))

    return mesh

def vol2voxelfe(voldata, templatefile, fileout, matprop=None, keywords=['NSET', 'ELSET'], voxelsize=[1, 1, 1], eltype='C3D8', matpropbits=8, verbose=False):
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
                  - NODES_S: Nodes on SOUTH surface of 3D model.
                  - NODES_N: Nodes on NORTH surface of 3D model.
                  - NODES_E: Nodes on EAST surface of 3D model.
                  - NODES_W: Nodes on WEST surface of 3D model.
                  - NODES_T: Nodes on TOP surface of 3D model.
                  - NODES_B: Nodes on BOTTOM surface of 3D model.
                These node sets are available for boundary conditions definition.
        * 'ELSET':
                Create boundary element sets. (Default = ON)
                If 'ELSET' is specified, the following element sets are created:
                  - ELEMS_S: Elements of SOUTH surface of 3D model.
                  - ELEMS_N: Elements of NORTH surface of 3D model.
                  - ELEMS_E: Elements of EAST surface of 3D model.
                  - ELEMS_W: Elements of WEST surface of 3D model.
                  - ELEMS_T: Elements of TOP surface of 3D model.
                  - ELEMS_B: Elements of BOTTOM surface of 3D model.
        * 'PROPERTY':
                Define an external material mapping law from template file. (Default = None)
                Use in combination with 'matprop' dictionary of material property files and corresponding GV ranges for the material mapping.
    voxelsize : float
        3D model voxelsize.
    eltype : str
        FE element type. The default is eight-node brick element (C3D8 and F3D8). See CalculiX node convention (Par. 6.2.1) at: http://www.dhondt.de/ccx_2.15.pdf
    matpropbits : int
        Bit depth for material mapping.
    """

    # verbose output
    if verbose:
        logging.basicConfig(level=logging.INFO)

    # get 3D dataset shape
    data_shape = voldata.shape

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
    logging.info('Detecting boundary nodes')
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
                    nodes_XYZ[node_i] = (voxelsize[0] * row, voxelsize[1] * col, voxelsize[2] * slice)

    # write ABAQUS *.inp output file ############################################
    logging.info('Start writing INP file')
    # open ABAQUS *.INP output file
    INP = open(fileout, 'w')

    # HEADER:
    INP.write('** ---------------------------------------------------------\n')
    INP.write('** Abaqus .INP file written by Voxel_FEM script on {}\n'.format(datetime.now()))
    INP.write('** ---------------------------------------------------------\n')
    # INP.write('*HEADING\n')
    # INP.write('main input {0}\n'.format(filename_in))

    # NODAL COORDINATES:
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
                    INP.write(
                        '{0:10d}, {n[0]:12.6f}, {n[1]:12.6f}, {n[2]:12.6f}\n'.format(node_i, n=nodes_XYZ[node_i]))

    # ELEMENTS AND ELEMENT SETS:
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

    # NODE SETS:
    if 'NSET' in keywords:
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
        logging.info('Additional node sets generated: NODES_S, NODES_N, NODES_E, NODES_W, NODES_T, NODES_B')

    # ELEMENT SETS:
    if 'ELSET' in keywords:
        INP.write('** Additional elset from voxel model. New generated elsets are:\n')
        INP.write('** ELEMS_S, ELEMS_N, ELEMS_E, ELEMS_W, ELEMS_T, ELEMS_B\n')

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
        logging.info('Additional element sets generated: ELEMS_S, ELEMS_N, ELEMS_E, ELEMS_W, ELEMS_T, ELEMS_B')

    # MATERIAL MAPPING:
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

    # copy line by line info on model solution and boundary conditions from Abaqus analysis template file
    try:
        template = open(templatefile, 'r')
    except IOError('Abaqus template file {} not found.'.format(templatefile)):
        exit(1)
    logging.info('Reading Abaqus template file {}'.format(templatefile))

    for line in template.readlines():
        # copy line to output Abaqus file
        INP.write('{}'.format(line))

    template.close()
    INP.close()
    logging.info(
        'Model with {0} nodes and {1} elements written to file {fname}'.format(len(nodes), n_els, fname=fileout))

def mesh2tetrafe(meshdata, templatefile, fileout, keywords=['NSET', 'ELSET'], float_fmt='.6e', verbose=False):
    """Generate ABAQUS tetrahedra Finite Element (FE) input file from 3D mesh. The output can be solved using ABAQUS or CalculiX.

    Boundary conditions (BCs), simulation steps and associated output requests are defined in a separate template file.
    See the templates contained in the folder "input_templates" for examples.

    Parameters
    ----------
    meshdata : meshio

    """

    # verbose output
    if verbose is True:
        logging.basicConfig(level=logging.INFO)

    # load 3D mesh using meshio module
    logging.info('Converting mesh with fields:\n')
    logging.info(vars(meshdata))

    if 'NSET' in keywords:
        # find model boundaries
        model_coors_max = np.amax(meshdata.points, 0)
        model_coors_min = np.amin(meshdata.points, 0)

        # volume extent in x, y and z
        extent = model_coors_max-model_coors_min

        # add dictionary of boundary point sets
        meshdata.point_sets = {
            'NODES_N': np.where(meshdata.points[:, 1] == model_coors_max[1])[0],
            'NODES_S': np.where(meshdata.points[:, 1] == model_coors_min[1])[0],
            'NODES_E': np.where(meshdata.points[:, 0] == model_coors_max[0])[0],
            'NODES_W': np.where(meshdata.points[:, 0] == model_coors_min[0])[0],
            'NODES_T': np.where(meshdata.points[:, 2] >= (model_coors_max[2]-0.01*abs(extent[2])))[0],
            'NODES_B': np.where(meshdata.points[:, 2] <= (model_coors_min[2]+0.01*abs(extent[2])))[0]
        }

    if 'ELSET' in keywords:
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

        for node_e in meshdata.point_sets['NODES_E']:
            ELEMS_E = np.append(ELEMS_E, np.where(np.any(meshdata.cells[0][1] == node_e, axis=1)))

        for node_w in meshdata.point_sets['NODES_W']:
            ELEMS_W = np.append(ELEMS_W, np.where(np.any(meshdata.cells[0][1] == node_w, axis=1)))

        for node_s in meshdata.point_sets['NODES_S']:
            ELEMS_S = np.append(ELEMS_S, np.where(np.any(meshdata.cells[0][1] == node_s, axis=1)))

        for node_n in meshdata.point_sets['NODES_N']:
            ELEMS_N = np.append(ELEMS_N, np.where(np.any(meshdata.cells[0][1] == node_n, axis=1)))

        for node_t in meshdata.point_sets['NODES_T']:
            ELEMS_T = np.append(ELEMS_T, np.where(np.any(meshdata.cells[0][1] == node_t, axis=1)))

        for node_b in meshdata.point_sets['NODES_B']:
            ELEMS_B = np.append(ELEMS_B, np.where(np.any(meshdata.cells[0][1] == node_b, axis=1)))

        meshdata.cell_sets = {'SET1': [np.arange(0, len(meshdata.cells[0][1]))]}

        if not ELEMS_N.size == 0:
            meshdata.cell_sets['ELEMS_N'] = [np.unique(ELEMS_N)]

        if not ELEMS_S.size == 0:
            meshdata.cell_sets['ELEMS_S'] = [np.unique(ELEMS_S)]

        if not ELEMS_W.size == 0:
            meshdata.cell_sets['ELEMS_W'] = [np.unique(ELEMS_W)]

        if not ELEMS_E.size == 0:
            meshdata.cell_sets['ELEMS_E'] = [np.unique(ELEMS_E)]

        if not ELEMS_T.size == 0:
            meshdata.cell_sets['ELEMS_T'] = [np.unique(ELEMS_T)]

        if not ELEMS_B.size == 0:
            meshdata.cell_sets['ELEMS_B'] = [np.unique(ELEMS_B)]

    else:
        meshdata.cell_sets = {'SET1': [np.arange(0, len(meshdata.cells[0][1]))]}

    # write Abaqus mesh using meshio
    meshio.abaqus.write(fileout, meshdata, float_fmt)

    # Add analysis definition to the end of the Abaqus INP from template
    # Open ABAQUS *.INP output file in Append mode
    INP = open(fileout, 'a')

    # copy line by line info on model solution and boundary conditions from Abaqus template file
    try:
        BCfile = open(templatefile, 'r')
    except IOError('Abaqus template file {} not found.'.format(templatefile)):
        exit(1)
    logging.info('Reading Abaqus template file {}'.format(templatefile))

    for line in BCfile.readlines():
        # copy line to output Abaqus file
        INP.write('{}'.format(line))

    BCfile.close()
    INP.close()
    logging.info('Data written to file {}'.format(fileout))

def matpropdictionary(proplist):
    """Compose dictionary of material properties and property mapping GV ranges.

    Parameters
    ----------
    proplist
        List of material property files followed by the corresponding GC range for material mapping.

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
            raise IOError('Incorrect use of --mapping flag: each material property file must be followed by MIN and MAX of the corresponding GV range.')

        else:
            matprop = {
                'file': [],
                'range': []
            }
            for prop in range(0, int(len(proplist)/3)):
                matprop['file'].append(proplist[3*prop])
                matprop['range'].append([int(proplist[3*prop+1]), int(proplist[3*prop+2])])

            return matprop

def main():
    description = textwrap.dedent('''\
                Script to generate ABAQUS Finite Element (FE) input files from 3D voxel data.

                The output of this script is an input file (.INP) in ABAQUS syntax
                that can be solved using ABAQUS or CALCULIX.

                The input 3D image is segmented and processed to maintain only the largest connected structure.
                3D voxel or tetrahedra FE input files can be generated.
                Meshes for 3D rendering of the outer model shell and of the whole selected volume can be generated.

                Only for voxel FE, the user can define a material mapping strategy
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

    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filein', type=str, help='<Required> Input filename (voxel data).')
    parser.add_argument('fileout', type=str, default=None, help='<Required> Output filename (Abaqus .INP).')
    parser.add_argument('-vs', '--voxelsize', type=float, default=[1., 1., 1.], nargs='+', help='Voxel size.')
    parser.add_argument('-r', '--resampling', type=float, default=1., help='Resampling factor.')
    parser.add_argument('-t', '--threshold', type=int, default=None, help='Threshold value.')
    parser.add_argument('--smooth', dest='smooth', action='store_true', help='Smooth image with gaussian filter before thresholding.')
    parser.add_argument('--caps', type=int, default=None, help='Add caps of given thickness to the bottom and top of the model for mesh creation.')
    parser.add_argument('--caps_val', type=int, default=0, help='Caps grey value.')
    parser.add_argument('--shell_mesh', dest='shell_mesh', action='store_true', help='Write VTK mesh of outer shell generated with PyMCubes.')
    parser.add_argument('--vol_mesh', dest='vol_mesh', action='store_true', help='Write VTK volume mesh of tetrahedra with pygalmesh.')
    parser.add_argument('--max_facet_distance', type=float, default=None, help='CGAL mesh parameter.')
    parser.add_argument('--max_cell_circumradius', type=float, default=None, help='CGAL mesh parameter.')
    parser.add_argument('--voxelfe', dest='voxelfe', action='store_true', help='Write voxel FE model (.INP) file.')
    parser.add_argument('--template', type=str, default=None, help='<Required by --voxelfe> Abaqus analysis template file (.INP).')
    parser.add_argument('-m', '--mapping', default=None, nargs='+', help='Template file for material property mapping. If more than one property is given, each property filename must followed by the corresponding GV range.')
    parser.add_argument('--tetrafe', dest='tetrafe', action='store_true', help='Write linear tetrahedra FE model (.INP) file.')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose output.')
    parser.set_defaults(smooth=False, shell_mesh=False, vol_mesh=False, voxelfe=False, tetrafe=False, verbose=False)

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    # filename base
    [fileout_base, type] = os.path.splitext(args.fileout)

    # isotropic voxel size (spacing) [mm]
    sp = np.array(args.voxelsize)

    # Read tiff stack #######################################################
    I = ru.read_tiff_stack(args.filein)

    # Gaussian smooth #######################################################
    if args.smooth:
        logging.info('Smoothing image with gaussian kernel')
        I = gaussian(I, sigma=1, preserve_range=True)

    # Resize the dataset ####################################################
    # use scikit
    if not args.resampling == 1:
        logging.info("Resampling input dataset with factor 1/{}".format(args.resampling))

        # resize the 3D data using spline interpolation of order 2
        I = ndimage.zoom(I, 1/args.resampling, output=None, order=2)

        # correct voxelsize
        sp = sp * args.resampling
        logging.info("New voxelsize: {}".format(sp))

    if args.verbose:
        # plot the midplanes through the stack
        ru.plot_midplanes(I)
        # write midplanes as .PNG
        ru.writemidplanes(I, fileout_base+".png")

    # Add image caps before thresholding ####################################
    if args.caps is not None:
        # add cap so that the mesh is closed
        I = ru.add_cap(I, cap_thickness=args.caps, cap_val=args.caps_val)

    # Binarise the dataset ##################################################
    # we do want bone = 1 and background = 0;
    if args.threshold is None:
        # use Otsu if threshold input not specified
        T = threshold_otsu(I)
        logging.info("Threshold value not given. Using Otsu method..")
        logging.info("Threshold: {}".format(T))

        if args.verbose:
            # plot the input image histogram
            fig2, ax2 = plt.subplots()
            plt.hist(I.ravel(), bins=100)
            plt.show()
            fig2.savefig(fileout_base+"_hist.png")
    else:
        T = int(args.threshold)
        logging.info("Threshold: {}".format(T))

    # apply the threshold
    BW = I > T

    # Detect the isolated cluster of voxels #################################
    L = remove_unconnected(BW)

    # Visualization with Napari #############################################
    # import napari
    # viewer = napari.view_image(L)

    # Generate outer shell mesh #############################################
    if args.shell_mesh:
        logging.info("Writing triangle mesh of outer shell")
        # shell mesh using pymcubes (high resolution mesh; caps excluded)
        vertices, triangles, shellmesh = shell_mesh(L, method='pymcubes')

        # write VTK mesh with meshio
        meshio.write_points_cells(fileout_base+"_shell.vtk", vertices.tolist(), [("triangle", triangles.tolist())])

        # shell mesh using pygalmesh (CGAL) (low resolution mesh; caps included)
        # check CGAL parameters
        # max_facet_distance, max_cell_circumradius = check_cgal_params(args.max_facet_distance, args.max_cell_circumradius, sp)

        # vertices, triangles, shellmesh = shell_mesh(L, 'pygalmesh', sp, max_facet_distance, max_cell_circumradius)
        # shellmesh.write(fileout_base + "_shell.vtk")

    # Generate volume mesh ##################################################
    if args.vol_mesh:
        # check CGAL parameters
        max_facet_distance, max_cell_circumradius = check_cgal_params(args.max_facet_distance, args.max_cell_circumradius, sp)

        # generate mesh with CGAL
        volmesh = cgal_mesh(L, sp, 'tetra', max_facet_distance, max_cell_circumradius)

        # write the mesh to file
        logging.info("Writing tetrahedra volume mesh")
        volmesh.write(fileout_base + "_vol.vtk")

    # Generate voxel FE model ##################################################
    if args.voxelfe:
        # check voxelfe args
        if not args.template:
            raise IOError('Analysis template file (--template) required with --voxelfe.')

        # build dictionary of material properties
        matprop = matpropdictionary(args.mapping)

        if args.mapping is None:
            # voxelFE of binary volume data; the material property definition is assumed to be in the analysis template file
            vol2voxelfe(L, args.template, fileout_base + "_voxelFE.inp", keywords=['NSET', 'ELSET'], voxelsize=sp, verbose=args.verbose)

        else:
            # mask the data (keep only the largest connected volume)
            masked = I
            masked[~L] = 0
            # voxelFE of greyscale volume data; material mapping on
            vol2voxelfe(masked, args.template, fileout_base + "_voxelFE.inp", matprop, keywords=['NSET', 'ELSET', 'PROPERTY'], voxelsize=sp, verbose=args.verbose)

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
                max_facet_distance, max_cell_circumradius = check_cgal_params(args.max_facet_distance, args.max_cell_circumradius, sp)

                # generate volume mesh
                volmesh = cgal_mesh(L, sp, 'tetra', max_facet_distance, max_cell_circumradius)

            elif volmesh is None:
                # check CGAL parameters
                max_facet_distance, max_cell_circumradius = check_cgal_params(args.max_facet_distance, args.max_cell_circumradius, sp)

                # generate volume mesh
                volmesh = cgal_mesh(L, sp, 'tetra', max_facet_distance, max_cell_circumradius)

            # tetraFE of the already meshed volume; the material property definition is assumed to be in the analysis template file
            mesh2tetrafe(volmesh, args.template, fileout_base + "_tetraFE.inp", keywords=['NSET', 'ELSET'], verbose=args.verbose)

        else:
            logging.exception("--tetrafe with mapping option not implemented yet.")

    return

if __name__ == '__main__':
    main()


