#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ciclope module for tetrahedra Finite Element model generation
"""

import logging
import numpy as np
import meshio
import mcubes
import pygalmesh
import math

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

def find_unique_edges(mesh):
    """
    Find unique edges in a tetrahedral mesh.

    This function iterates through all the tetrahedral cells in the mesh
    and extracts unique edges by considering each pair of vertices in the tetrahedra.
    Each edge is represented as a tuple of two sorted vertex indices.

    Parameters
    ----------
    mesh : Mesh
        A mesh object that contains cells, where each cell represents a tetrahedron.

    Returns
    -------
    unique_edges : set of tuple of int
        A set containing unique edges, where each edge is represented as a
        tuple of two vertex indices (sorted in ascending order). Each tuple
        represents an undirected edge between two vertices in the mesh.
    """

    unique_edges = set()
    for cell in mesh.cells:
        if cell.type == "tetra":
            for tetra in cell.data:
                for i in range(4):
                    for j in range(i+1, 4):
                        edge = tuple(sorted([tetra[i], tetra[j]]))
                        unique_edges.add(edge)

    return unique_edges

def add_midpoints_to_mesh(mesh, sorted_edges):
    """
    Add midpoints of edges to the mesh and update edge-to-midpoint mapping.

    This function computes the midpoints of the given edges, adds these midpoints
    to the mesh's points, and updates a dictionary that maps each edge to its 
    corresponding midpoint index in the mesh.

    Parameters
    ----------
    mesh : Mesh
        The mesh object containing points and cells.
    sorted_edges : list of tuple of int
        A list of sorted edges, where each edge is represented as a tuple 
        containing two vertex indices (sorted in ascending order).

    Returns
    -------
    edge_to_midpoint : dict of tuple of int to int
        A dictionary mapping each edge (represented as a tuple of two vertex 
        indices) to the index of its midpoint in the mesh's points array.
    """

    edge_num = len(sorted_edges)
    edge_to_midpoint = {}

    # Preallocate a 2D array for the new midpoints
    new_points = np.empty((edge_num, 3))
    old_point_num = len(mesh.points)

    for i, edge in enumerate(sorted_edges):
        node1, node2 = edge
        # Calculate the midpoint and store it in the new_points array
        new_points[i] = (mesh.points[node1] + mesh.points[node2]) / 2
        # Map the edge to the index of the new midpoint
        edge_to_midpoint[edge] = old_point_num + i

    # Concatenate the new points with the existing points in the mesh
    mesh.points = np.concatenate((mesh.points, new_points))

    return edge_to_midpoint

def create_new_cells(mesh, edge_to_midpoint):
    """
    Create a new list of tetrahedral cells with midpoint nodes.

    This function generates new cells for the mesh by adding midpoint nodes
    to the edges of existing tetrahedral cells. Each new cell includes the
    original vertices and the additional midpoint nodes, resulting in a 
    10-node tetrahedral element (tetra10).

    Parameters
    ----------
    mesh : Mesh
        The mesh object containing cells and points.
    edge_to_midpoint : dict of tuple of int to int
        A dictionary mapping each edge (represented as a tuple of two vertex
        indices) to the index of its corresponding midpoint in the mesh's points
        array.

    Returns
    -------
    new_cells : list of meshio.CellBlock
        A list containing the new cells with the added midpoint nodes. Each cell
        is represented as a `meshio.CellBlock` object with 'tetra10' type, and the
        `data` attribute contains a list of vertex indices, including the midpoints.
    """

    new_cells_data = []

    for i, cell in enumerate(mesh.cells):
        if cell.type == 'tetra':
            for nodes in cell.data:
                new_cell = list(nodes)

                edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
                for edge in edges:
                    key = tuple(sorted((nodes[edge[0]], nodes[edge[1]])))
                    new_cell.append(edge_to_midpoint[key])

                new_cells_data.append(new_cell)

    new_cells = [meshio.CellBlock('tetra10', np.array(new_cells_data))]
    mesh.cells = new_cells

    return new_cells

def mesh2tetrafe(meshdata, templatefile, fileout, keywords=['NSET', 'ELSET'], float_fmt='.6e', verbose=False, bound_tol=None):
    """
    Generate an ABAQUS tetrahedral Finite Element (FE) input file from a 3D mesh.

    This function converts a mesh (typically created with meshio) into an ABAQUS input file (.inp),
    which can be used for simulations with ABAQUS or CalculiX. It creates boundary node sets (NSET)
    and, if specified, boundary element sets (ELSET) based on geometric criteria. Additionally, the function
    appends analysis definitions and boundary condition specifications to the end of the file by copying them
    from a template file.

    Boundary conditions, simulation steps, and associated output requests are defined in a separate template file.
    See the templates in the "input_templates" folder for practical examples.

    Parameters
    ----------
    meshdata : meshio.Mesh
        Mesh data including node coordinates, element connectivity, and optionally cell data.
        This object is typically generated by reading a mesh file with meshio.
    templatefile : str
        Path to the ABAQUS template file containing the analysis definition (e.g., simulation steps, boundary conditions).
    fileout : str
        Name of the output file (.inp) where the generated ABAQUS input file will be written.
    keywords : list of str, optional
        List of keywords specifying which boundary sets to create. Supported keywords are:
          - 'NSET': Create boundary node sets. When specified, the following node sets are created:
                * NODES_X0: Nodes on the West (X-) side of the model.
                * NODES_X1: Nodes on the East (X+) side of the model.
                * NODES_Y0: Nodes on the South (Y-) side of the model.
                * NODES_Y1: Nodes on the North (Y+) side of the model.
                * NODES_Z0: Nodes on the Bottom (Z-) side of the model.
                * NODES_Z1: Nodes on the Top (Z+) side of the model.
                Additionally, example corner node sets (e.g., NODES_X0Y0Z0 and NODES_X0Y0Z1) are created.
          - 'ELSET': Create boundary element sets. When specified, for each boundary node set a corresponding
                element set is created that includes all elements containing at least one node from the respective set.
                Moreover, if the meshdata contains cell data under the key 'medit:ref', additional element sets are
                generated based on the unique values found.
    float_fmt : str, optional
        Format string for floating-point numbers when writing the ABAQUS file.
        The default value is '.6e'.
    verbose : bool, optional
        If True, enables detailed logging output during execution, which is useful for debugging and tracking the process.
    bound_tol : float, list of float, or np.array, optional
        Tolerance used to determine if a node lies on the boundary of the model.
        - If None (default), it is set to 1% of the model's extent in each direction.
        - If a single float is provided, the same tolerance is applied to all three directions (X, Y, Z).
        - If a list or array is provided, it must contain 1 or 3 values corresponding to the tolerances for the X, Y, and Z directions.

    Returns
    -------
    None
        The function writes the generated ABAQUS input file to the path specified by 'fileout'
        and does not return any value.

    Notes
    -----
    The generated ABAQUS input file is composed of two parts:
      1. The mesh section, which includes the node coordinates and element connectivity formatted according
         to ABAQUS specifications using meshio.
      2. The analysis definition section, which is appended to the end of the file by copying the contents
         from the provided template file.
    """

    if verbose:
        logging.basicConfig(level=logging.INFO)

    # Calcola i limiti del modello
    model_coors_max = np.amax(meshdata.points, axis=0)
    model_coors_min = np.amin(meshdata.points, axis=0)
    extent = model_coors_max - model_coors_min

    # Tolleranza di confine: se non specificata, 1% dell'estensione del modello
    if bound_tol is None:
        bound_tol = 0.01 * extent
    else:
        bound_tol = np.array(bound_tol)
        if bound_tol.size == 1:
            bound_tol = np.repeat(bound_tol, 3)
        elif bound_tol.size != 3:
            raise ValueError("bound_tol deve essere un float, una lista o un array con 1 o 3 valori")

    # Definizione dei set di nodi di confine (NSET)
    if 'NSET' in keywords:
        meshdata.point_sets = {
            'NODES_X0': np.where(meshdata.points[:, 0] <= (model_coors_min[0] + bound_tol[0]))[0],            
            'NODES_X1': np.where(meshdata.points[:, 0] >= (model_coors_max[0] - bound_tol[0]))[0],
            'NODES_Y0': np.where(meshdata.points[:, 1] <= (model_coors_min[1] + bound_tol[1]))[0],
            'NODES_Y1': np.where(meshdata.points[:, 1] >= (model_coors_max[1] - bound_tol[1]))[0],
            'NODES_Z0': np.where(meshdata.points[:, 2] <= (model_coors_min[2] + bound_tol[2]))[0],
            'NODES_Z1': np.where(meshdata.points[:, 2] >= (model_coors_max[2] - bound_tol[2]))[0]
        }

    # Aggiunge i set per i nodi degli angoli (esempio)
    meshdata.point_sets['NODES_X0Y0Z0'] = meshdata.point_sets['NODES_Z0'][0:2]
    meshdata.point_sets['NODES_X0Y0Z1'] = meshdata.point_sets['NODES_Z1'][0:2]

    # Definizione dei set di elementi (ELSET)
    if 'ELSET' in keywords:
        # Inizializza il set che include tutti gli elementi
        cell_sets = {'SETALL': [np.arange(len(meshdata.cells[0][1]))]}
        
        # Mappa tra set di elementi e set di nodi di confine corrispondenti
        mapping = {
            'ELEMS_X0': 'NODES_X0',
            'ELEMS_X1': 'NODES_X1',
            'ELEMS_Y0': 'NODES_Y0',
            'ELEMS_Y1': 'NODES_Y1',
            'ELEMS_Z0': 'NODES_Z0',
            'ELEMS_Z1': 'NODES_Z1'
        }
        
        # Array contenente le celle (ogni riga sono gli indici dei nodi che compongono l'elemento)
        cells_array = meshdata.cells[0][1]
        
        # Se esistono dati per le celle, crea ulteriori set in base ai valori unici
        if hasattr(meshdata, 'cell_data'):
            cell_data_unique = np.unique(meshdata.cell_data['medit:ref'])
            for val in cell_data_unique:
                indices = np.where(meshdata.cell_data['medit:ref'] == val)[0]
                cell_sets['SET' + str(val)] = [indices]

        # Per ogni set di nodi, utilizza la vectorizzazione per trovare rapidamente gli elementi
        for cell_set_name, node_set_name in mapping.items():
            node_set = meshdata.point_sets[node_set_name]
            # Crea una maschera: per ogni elemento controlla se contiene almeno un nodo del set
            mask = np.any(np.isin(cells_array, node_set), axis=1)
            elem_indices = np.where(mask)[0]
            if elem_indices.size > 0:
                # Salva gli indici unici degli elementi trovati
                cell_sets[cell_set_name] = [np.unique(elem_indices)]
        
        meshdata.cell_sets = cell_sets

    else:
        meshdata.cell_sets = {'SETALL': [np.arange(len(meshdata.cells[0][1]))]}

    # Scrive la mesh in formato ABAQUS usando meshio
    meshio.abaqus.write(fileout, meshdata, float_fmt)

    # Aggiunge la definizione dell'analisi alla fine del file di input ABAQUS, copiandolo da un template
    try:
        with open(templatefile, 'r') as BCfile, open(fileout, 'a') as INP:
            for line in BCfile:
                INP.write(line)
    except IOError:
        logging.error("Abaqus template file {} non trovato.".format(templatefile))
        exit(1)