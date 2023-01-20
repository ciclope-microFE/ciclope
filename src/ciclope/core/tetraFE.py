#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ciclope module for tetrahedra Finite Element model generation
"""

import logging
import numpy as np
import meshio
import mcubes

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
        Mesh data.
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
    mesh : meshio
        Mesh data.

    """

    import pygalmesh
    # generate mesh
    mesh = pygalmesh.generate_from_array(np.transpose(bwimage, [2, 1, 0]).astype('uint8'), tuple(voxelsize), max_facet_distance=max_facet_distance, max_cell_circumradius=max_cell_circumradius, verbose=False)

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

def mesh2tetrafe(meshdata, templatefile, fileout, keywords=['NSET', 'ELSET'], float_fmt='.6e', verbose=False):
    """Generate ABAQUS tetrahedra Finite Element (FE) input file from 3D mesh. The output can be solved using ABAQUS or CalculiX.

    Boundary conditions (BCs), simulation steps and associated output requests are defined in a separate template file.
    See the templates contained in the folder "input_templates" for examples.

    Parameters
    ----------
    meshdata : meshio
        Mesh data.
    templatefile : str
        Analysis template file.
    fileout : str
        Output .INP file.
    keywords : str
        SUPPORTED ABAQUS KEYWORDS:

        For a list of all Abaqus keywords and their description visit:
        https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-c-gen-kwbrowser.htm#simacae-c-gen-kwbrowser__simacae-gen-xsl-U

        * 'NSET':
                Create boundary node sets. (Default = ON)
                If 'NSET' is specified, the following node sets are created:
                  - NODES_X0: Nodes on WEST (X-) surface of 3D model.
                  - NODES_X1: Nodes on EAST (X+) surface of 3D model.
                  - NODES_Y0: Nodes on SOUTH (Y-) surface of 3D model.
                  - NODES_Y1: Nodes on NORTH (Y+) surface of 3D model.
                  - NODES_Z0: Nodes on BOTTOM (Z-) surface of 3D model.
                  - NODES_Z1: Nodes on TOP (Z+) surface of 3D model.
                  - NODES_X0Y0Z0: 2 nodes on (0,0,0) model corner.
                  - NODES_X0Y0Z1: 2 nodes on (0,0,1) model corner.
                These node sets are available for boundary conditions definition.
        * 'ELSET':
                Create boundary element sets. (Default = ON)
                If 'ELSET' is specified, the following element sets are created:
                  - ELEMS_X0: Elements of WEST (X-) surface of 3D model.
                  - ELEMS_X1: Elements of EAST (X+) surface of 3D model.
                  - ELEMS_Y0: Elements of SOUTH (Y-) surface of 3D model.
                  - ELEMS_Y1: Elements of NORTH (Y+) surface of 3D model.
                  - ELEMS_Z0: Elements of BOTTOM (Z-) surface of 3D model.
                  - ELEMS_Z1: Elements of TOP (Z+) surface of 3D model.
    float_fmt : float
        Precision for Abaqus input file writing.
    verbose : bool
        Verbose output.
    """

    # verbose output
    if verbose is True:
        logging.basicConfig(level=logging.INFO)

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
            'NODES_Y1': np.where(meshdata.points[:, 1] == model_coors_max[1])[0],
            'NODES_Y0': np.where(meshdata.points[:, 1] == model_coors_min[1])[0],
            'NODES_X1': np.where(meshdata.points[:, 0] == model_coors_max[0])[0],
            'NODES_X0': np.where(meshdata.points[:, 0] == model_coors_min[0])[0],
            'NODES_Z1': np.where(meshdata.points[:, 2] >= (model_coors_max[2]-0.01*abs(extent[2])))[0],
            'NODES_Z0': np.where(meshdata.points[:, 2] <= (model_coors_min[2]+0.01*abs(extent[2])))[0]
        }

    meshdata.point_sets['NODES_X0Y0Z0'] = meshdata.point_sets['NODES_Z0'][0:2]
    meshdata.point_sets['NODES_X0Y0Z1'] = meshdata.point_sets['NODES_Z1'][0:2]

    if 'ELSET' in keywords:
        # dictionary of boundary element sets
        elset = {
            'ELEMS_Y1': [],
            'ELEMS_Y0': [],
            'ELEMS_X0': [],
            'ELEMS_X1': [],
            'ELEMS_Z1': [],
            'ELEMS_Z0': []
        }

        # find boundary cells
        ELEMS_X1 = np.array([]).astype('int')
        ELEMS_X0 = np.array([]).astype('int')
        ELEMS_Y0 = np.array([]).astype('int')
        ELEMS_Y1 = np.array([]).astype('int')
        ELEMS_Z1 = np.array([]).astype('int')
        ELEMS_Z0 = np.array([]).astype('int')

        for node_e in meshdata.point_sets['NODES_X1']:
            ELEMS_X1 = np.append(ELEMS_X1, np.where(np.any(meshdata.cells[0][1] == node_e, axis=1)))

        for node_w in meshdata.point_sets['NODES_X0']:
            ELEMS_X0 = np.append(ELEMS_X0, np.where(np.any(meshdata.cells[0][1] == node_w, axis=1)))

        for node_s in meshdata.point_sets['NODES_Y0']:
            ELEMS_Y0 = np.append(ELEMS_Y0, np.where(np.any(meshdata.cells[0][1] == node_s, axis=1)))

        for node_n in meshdata.point_sets['NODES_Y1']:
            ELEMS_Y1 = np.append(ELEMS_Y1, np.where(np.any(meshdata.cells[0][1] == node_n, axis=1)))

        for node_t in meshdata.point_sets['NODES_Z1']:
            ELEMS_Z1 = np.append(ELEMS_Z1, np.where(np.any(meshdata.cells[0][1] == node_t, axis=1)))

        for node_b in meshdata.point_sets['NODES_Z0']:
            ELEMS_Z0 = np.append(ELEMS_Z0, np.where(np.any(meshdata.cells[0][1] == node_b, axis=1)))

        meshdata.cell_sets = {'SETALL': [np.arange(0, len(meshdata.cells[0][1]))]}

        # if cell data exists create element sets based on it
        if hasattr(meshdata, 'cell_data'):

            # find unique cell data
            cell_data_unique = np.unique(meshdata.cell_data['medit:ref'])

            # for each unique cell scalar find indexes of the corresponding cells and create element set of them
            for val in cell_data_unique:
                meshdata.cell_sets['SET' + str(val)] = [np.where(meshdata.cell_data['medit:ref'] == val)[0]]
                # print(val)

        if not ELEMS_Y1.size == 0:
            meshdata.cell_sets['ELEMS_Y1'] = [np.unique(ELEMS_Y1)]

        if not ELEMS_Y0.size == 0:
            meshdata.cell_sets['ELEMS_Y0'] = [np.unique(ELEMS_Y0)]

        if not ELEMS_X0.size == 0:
            meshdata.cell_sets['ELEMS_X0'] = [np.unique(ELEMS_X0)]

        if not ELEMS_X1.size == 0:
            meshdata.cell_sets['ELEMS_X1'] = [np.unique(ELEMS_X1)]

        if not ELEMS_Z1.size == 0:
            meshdata.cell_sets['ELEMS_Z1'] = [np.unique(ELEMS_Z1)]

        if not ELEMS_Z0.size == 0:
            meshdata.cell_sets['ELEMS_Z0'] = [np.unique(ELEMS_Z0)]

    else:
        meshdata.cell_sets = {'SETALL': [np.arange(0, len(meshdata.cells[0][1]))]}

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

