#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ciclope module for voxel Finite Element model generation
"""

import os
import h5py
import logging
import numpy as np
from datetime import datetime
from ciclope.utils.preprocess import _collect_from_elements
from ciclope.utils.recon_utils import bbox

def vol2ugrid(voldata, voxelsize=[1, 1, 1], GVmin=0, plane_lock_num=1, node_level_lock_num=1, locking_strategy="plane", refnodes=False, verbose=False):
    """
    Generate unstructured grid mesh from 3D volume data.

    Parameters
    ----------
    voldata : ndarray
        3D voxel data.
    voxelsize : list or array
        3D model voxelsize.
    GVmin : int or float
        Minimum Grey Value considered for meshing. By default, GVmin=0: all zeroes are considered background.
    plane_lock_num : int
        Number of volume slices to consider for locking in "plane" mode.
    node_level_lock_num : int
        Number of node levels to consider for locking in "exact" mode.
    locking_strategy : str
        Strategy for node locking:
            - "plane" locks the entire strip of cells (default), 
        - "exact" locks nodes exactly at specific levels.
    refnodes : bool
        Return dictionary of reference nodes on the model boundaries. Even if this option is not activated, the returned
        mesh will contain the following nodes and elements sets:
            - NODES_X0: Nodes on WEST (X-) surface of 3D model.
            - NODES_X1: Nodes on EAST (X+) surface of 3D model.
            - NODES_Y0: Nodes on SOUTH (Y-) surface of 3D model.
            - NODES_Y1: Nodes on NORTH (Y+) surface of 3D model.
            - NODES_Z0: Nodes on BOTTOM (Z-) surface of 3D model.
            - NODES_Z1: Nodes on TOP (Z+) surface of 3D model.
            - NODES_X0Y0Z0: 2 nodes on (0,0,0) model corner.
            - NODES_X0Y0Z1: 2 nodes on (0,0,1) model corner.

            - ELEMS_X0: Elements of WEST (X-) surface of 3D model.
            - ELEMS_X1: Elements of EAST (X+) surface of 3D model.
            - ELEMS_Y0: Elements of SOUTH (Y-) surface of 3D model.
            - ELEMS_Y1: Elements of NORTH (Y+) surface of 3D model.
            - ELEMS_Z0: Elements of BOTTOM (Z-) surface of 3D model.
            - ELEMS_Z1: Elements of TOP (Z+) surface of 3D model.
    verbose : bool
        Activate verbose output.

    Returns
    -------
    mesh : meshio.Mesh
        Unstructured grid mesh.
    refnodes_dict : dict (optional)
        centroids on the model boundaries (X0, X1, Y0, Y1, Z0, Z1) if refnodes is True.
    """
    if verbose:
        logging.basicConfig(level=logging.INFO)

    # Get 3D dataset shape
    data_shape = voldata.shape

    # voxelsize list if only one value is given
    if not isinstance(voxelsize, (list, np.ndarray)):
        voxelsize = [voxelsize, voxelsize, voxelsize]
    if len(voxelsize) == 1:
        voxelsize = [voxelsize[0]] * 3

    # Initialization of lists and dictionaries
    nodes = []         # Node coordinates
    cells = []         # List of cells (elements)
    cell_GV = []       # List of grey values for each cell

    # Dictionary of boundary node sets
    nodes_set = {
        'NODES_X0': [],
        'NODES_X1': [],
        'NODES_Y0': [],
        'NODES_Y1': [],
        'NODES_Z0': [],
        'NODES_Z1': [],
        'NODES_X0Y0Z0': [],
        'NODES_X0Y0Z1': []
    }
    
    # Dictionary of boundary cell sets
    cells_set = {
        'CELLS_X0': [],
        'CELLS_X1': [],
        'CELLS_Y0': [],
        'CELLS_Y1': [],
        'CELLS_Z0': [],
        'CELLS_Z1': []
    }

    # Reference nodes coordinates (one refnode per boundary)
    refnode_X0 = np.zeros(3)
    refnode_X1 = np.zeros(3)
    refnode_Y0 = np.zeros(3)
    refnode_Y1 = np.zeros(3)
    refnode_Z0 = np.zeros(3)
    refnode_Z1 = np.zeros(3)
    
    # Binary mask for existing nodes
    if locking_strategy == "exact":
        existing_nodes = np.full((data_shape[0]+1)*(data_shape[1]+1)*(data_shape[2]+1), False)
    
    # Dictionary to store indices of nodes belonging to the surfaces
    elif locking_strategy == "plane":    
        planes_nodes = {
            'NODES_Z0': set(),
            'NODES_Z1': set(),
            'NODES_X0': set(),
            'NODES_X1': set(),
            'NODES_Y0': set(),
            'NODES_Y1': set()
        }
        
    else:
        raise ValueError("locking_strategy must be 'plane' or 'exact'")

    # Preliminary calculations for the number of nodes per row and per slice
    row_nodes = data_shape[2] + 1                # Nodes along Y (columns + 1)
    slice_nodes = (data_shape[2] + 1) * (data_shape[1] + 1)  # Nodes per slice
    cell_i = 0

    logging.info('Calculating cell array')
    # Loop to create cells (elements) from the voxel matrix
    for slice in range(data_shape[0]):
        for row in range(data_shape[1]):
            for col in range(data_shape[2]):

                # Current voxel value
                GV = voldata[slice, row, col]
                if GV > GVmin:
                    cell_i += 1
                    # Calculate the node indices for the element (cell)
                    # See eight-node brick cell (C3D8 and F3D8) node convention (Par. 6.2.1) at: http://www.dhondt.de/ccx_2.15.pdf
                    cell_nodes = [
                        slice_nodes * slice + row_nodes * row + col,
                        slice_nodes * slice + row_nodes * row + col + 1,
                        slice_nodes * slice + row_nodes * (row + 1) + col + 1,
                        slice_nodes * slice + row_nodes * (row + 1) + col,
                        slice_nodes * (slice + 1) + row_nodes * row + col,
                        slice_nodes * (slice + 1) + row_nodes * row + col + 1,
                        slice_nodes * (slice + 1) + row_nodes * (row + 1) + col + 1,
                        slice_nodes * (slice + 1) + row_nodes * (row + 1) + col
                    ]
                    
                    # append to cell list
                    cells.append(cell_nodes)
                    
                    # append to cell_GV list
                    cell_GV.append(GV)

                    # Handle locking along Z depending on the selected strategy
                    if locking_strategy == "plane":
                        if slice < plane_lock_num:
                            cells_set['CELLS_Z0'].append(cell_i)
                            for i in cell_nodes:
                                planes_nodes['NODES_Z0'].add(i)
                        if slice >= data_shape[0] - plane_lock_num:
                            cells_set['CELLS_Z1'].append(cell_i)
                            for i in cell_nodes:
                                planes_nodes['NODES_Z1'].add(i)
                        if col < plane_lock_num:
                            cells_set['CELLS_X0'].append(cell_i)
                            for i in cell_nodes:
                                planes_nodes['NODES_X0'].add(i)
                        if col >= data_shape[2] - plane_lock_num:
                            cells_set['CELLS_X1'].append(cell_i)
                            for i in cell_nodes:
                                planes_nodes['NODES_X1'].add(i)
                        if row < plane_lock_num:
                            cells_set['CELLS_Y0'].append(cell_i)
                            for i in cell_nodes:
                                planes_nodes['NODES_Y0'].add(i)
                        if row >= data_shape[1] - plane_lock_num:
                            cells_set['CELLS_Y1'].append(cell_i)
                            for i in cell_nodes:
                                planes_nodes['NODES_Y1'].add(i)

                    elif locking_strategy == "exact":
                        # Store if the cell nodes exist
                        for i in cell_nodes:
                            existing_nodes[i] = True
    logging.info('Detecting node coordinates and boundary nodes')
    # Loop to generate coordinates for all nodes and build the sets of boundary nodes
    node_i = 0
    total_slices = data_shape[0] + 1
    total_rows = data_shape[1] + 1
    total_cols = data_shape[2] + 1
    for slice in range(total_slices):
        for row in range(total_rows):
            for col in range(total_cols):
                node_coors = [voxelsize[0] * col, voxelsize[1] * row, voxelsize[2] * slice]
                nodes.append(node_coors)
   
                if locking_strategy == "plane":
                    if node_i in planes_nodes['NODES_Z0']:
                        nodes_set['NODES_Z0'].append(node_i)
                        refnode_Z0 += node_coors
                    if node_i in planes_nodes['NODES_Z1']:
                        nodes_set['NODES_Z1'].append(node_i)
                        refnode_Z1 += node_coors
                    if node_i in planes_nodes['NODES_X0']:
                        nodes_set['NODES_X0'].append(node_i)
                        refnode_X0 += node_coors
                    if node_i in planes_nodes['NODES_X1']:
                        nodes_set['NODES_X1'].append(node_i)
                        refnode_X1 += node_coors
                    if node_i in planes_nodes['NODES_Y0']:
                        nodes_set['NODES_Y0'].append(node_i)
                        refnode_Y0 += node_coors
                    if node_i in planes_nodes['NODES_Y1']:
                        nodes_set['NODES_Y1'].append(node_i)
                        refnode_Y1 += node_coors
                        
                elif locking_strategy == "exact":
                    if existing_nodes[node_i]:
                        # compose dictionary of boundary node sets
                        if slice < node_level_lock_num:
                            nodes_set['NODES_Z0'].append(node_i)
                            refnode_Z0 += node_coors
                        if slice > data_shape[0]-node_level_lock_num:
                            nodes_set['NODES_Z1'].append(node_i)
                            refnode_Z1 += node_coors
                        if col < node_level_lock_num:
                            nodes_set['NODES_X0'].append(node_i)
                            refnode_X0 += node_coors
                        if col > data_shape[2]-node_level_lock_num:
                            nodes_set['NODES_X1'].append(node_i)
                            refnode_X1 += node_coors
                        if row < node_level_lock_num:
                            nodes_set['NODES_Y0'].append(node_i)
                            refnode_Y0 += node_coors
                        if row > data_shape[1]-node_level_lock_num:
                            nodes_set['NODES_Y1'].append(node_i)
                            refnode_Y1 += node_coors
                            
                node_i += 1
    
    # For completeness, also define the corners (using the first two nodes found on surface Z)
    nodes_set['NODES_X0Y0Z0'] = nodes_set['NODES_Z0'][0:2]
    nodes_set['NODES_X0Y0Z1'] = nodes_set['NODES_Z1'][0:2]

    # Reference nodes are the barycenters of boundary node sets
    refnode_X0 /= len(nodes_set['NODES_X0']) if nodes_set['NODES_X0'] else 1
    refnode_X1 /= len(nodes_set['NODES_X1']) if nodes_set['NODES_X1'] else 1
    refnode_Y0 /= len(nodes_set['NODES_Y0']) if nodes_set['NODES_Y0'] else 1
    refnode_Y1 /= len(nodes_set['NODES_Y1']) if nodes_set['NODES_Y1'] else 1
    refnode_Z0 /= len(nodes_set['NODES_Z0']) if nodes_set['NODES_Z0'] else 1
    refnode_Z1 /= len(nodes_set['NODES_Z1']) if nodes_set['NODES_Z1'] else 1

    refnodes_dict = {
        "X0": refnode_X0,
        "X1": refnode_X1,
        "Y0": refnode_Y0,
        "Y1": refnode_Y1,
        "Z0": refnode_Z0,
        "Z1": refnode_Z1,
    }

    # Generate meshio object
    import meshio
    cells_dict = [("hexahedron", cells)]
    mesh = meshio.Mesh(nodes, cells_dict, cell_data={"GV": [cell_GV]},
                        point_sets=nodes_set, cell_sets=cells_set)

    logging.info('Generated the following mesh with {0} nodes and {1} elements:'
                 .format(len(mesh.points), len(mesh.cells[0])))
    logging.info(mesh)

    if refnodes:
        return mesh, refnodes_dict
    else:
        return mesh
        
def vol2h5ParOSol(voldata, fileout, topDisplacement, voxelsize=1, poisson_ratio=0.3, young_modulus=18e3, topHorizontaFixedlDisplacement=True, plane_lock_num = 1, verbose=False):
    """Generate ParOSol HDF5 (.h5) input file from 3D volume data.
    Before to generate ParOSol HDF5 file, the Bounding BOX (bbox class) limits the input binary image.
    Info on HDF5 file type for ParOSol solver at: https://github.com/reox/parosol-tu-wien/blob/master/doc/file_format.md

    Parameters
    ----------
    
    voldata : ndarray
        3D voxel data.
    fileout : str
        Output .h5 file.
    topDisplacement: float
        Vertical displacement imposed at the top. 
    voxelsize : float
        3D model voxelsize.
    poisson_ratio: float
        Poisson ratio.
    Young_modulus: float
        Young's modulus [MPa].
    topHorizontaFixedlDisplacement: bool 
        if True X and Y displacements fixed at the top; if False no displacements fixed at the top.
    plane_lock_num: int
        Number of planes where boundary conditions (BCs) are applied.
    plane_lock_num : int
        Number of planes where boundary conditions (BCs) are applied. 
        By default, BCs are applied only to the first layer nodes at the bases.        
    verbose : bool
        Activate verbose output.
    """    
    #definitions
    H5T_IEEE_F64LE='<f8'
    H5T_IEEE_F32LE='<f4'
    H5T_STD_U16LE='<u2'
   
    # verbose output
    if verbose:
        logging.basicConfig(level=logging.INFO)
        
    #file and main group(Image_Data) creation
    logging.info('Opening Output file')
    file=h5py.File(fileout, 'w')

    logging.info('creating h5 Image_Data Group')
    imgData=file.create_group("Image_Data")
   
    # get 3D input dataset shape
    data_shape = voldata.shape
    
    #using bbox utility to calculate output dataset shape (excluding empty planes)
    [origin, dims]=bbox(voldata)
    
    #use the same format as data_shape for origin and dims (reorder from Y,X,Z to Z,Y,X)
    origin=[origin[2],origin[0],origin[1]]
    dims=[dims[2],dims[0],dims[1]]
    
    # Compute end indices along each axis
    end = [origin[i] + dims[i] for i in range(3)]
    
    # Define slicing ranges for sub-block extraction
    slice_start, slice_end = origin[0], end[0] + 1
    row_start,   row_end   = origin[1], end[1] + 1
    col_start,   col_end   = origin[2], end[2] + 1

    # Extract the relevant sub-block of voxel data
    sub_block = voldata[slice_start:slice_end,
                        row_start:row_end,
                        col_start:col_end]
    # Scale voxel values by Young's modulus to get material stiffness
    array_data = (sub_block * young_modulus).astype(np.float32)

    # Store image data and material parameters in HDF5
    imgData.create_dataset("Image", data=array_data, dtype=H5T_IEEE_F32LE)
    
    #setting voxel size
    logging.info('Setting voxel size')
    imgData.create_dataset("Voxelsize", data=voxelsize, dtype=H5T_IEEE_F64LE)

    #setting poisson ratio
    logging.info('Setting Poisson ratio')
    imgData.create_dataset("Poisons_ratio", data=poisson_ratio, dtype=H5T_IEEE_F64LE)
    
    # Initialize sets to collect boundary nodes for bottom (lb) and top (tb)
    lbNodesSet = set()
    tbNodesSet = set()

    # PLANE: select voxels in first/last plane_lock_num layers
    lowBoundaryEls = []
    topBoundaryEls = []
    for i in range(plane_lock_num):
        # Identify non-zero voxels in bottom plane i
        low_mask = (sub_block[i] != 0)
        low2d = np.argwhere(low_mask)
        # Prepend the plane index to each (y,x) coordinate
        low2d = np.hstack((np.full((low2d.shape[0],1), i), low2d))
        lowBoundaryEls.append(low2d)

        # Identify non-zero voxels in top plane (from end)
        idx = sub_block.shape[0] - 1 - i
        top_mask = (sub_block[idx] != 0)
        top2d = np.argwhere(top_mask)
        top2d = np.hstack((np.full((top2d.shape[0],1), idx), top2d))
        topBoundaryEls.append(top2d)
    # Concatenate all boundary element arrays
    lowBoundaryEls = np.vstack(lowBoundaryEls)
    topBoundaryEls = np.vstack(topBoundaryEls)

    # Bottom: fix X,Y,Z; Top: conditionally fix X,Y based on flag
    _collect_from_elements(lowBoundaryEls, lbNodesSet, (0,1,2))
    if topHorizontaFixedlDisplacement:
        _collect_from_elements(topBoundaryEls, tbNodesSet, (0,1,2))
    else:
        _collect_from_elements(topBoundaryEls, tbNodesSet, (2,))

    # Merge bottom+top sets and prepare arrays for HDF5
    all_nodes = list(lbNodesSet) + list(tbNodesSet)
    fixDispCoord = np.array(all_nodes, dtype=np.uint16)  # (n_nodes,4)

    # Initialize displacement values (0 for all), then set top Z displacements
    values = np.zeros(len(all_nodes), dtype=np.float32)
    for i, (_, _, _, d) in enumerate(all_nodes[len(lbNodesSet):], start=len(lbNodesSet)):
        if d == 2:
            values[i] = topDisplacement

    # Write Fixed_Displacement datasets
    imgData.create_dataset(
        "Fixed_Displacement_Coordinates",
        data=fixDispCoord,
        dtype=H5T_STD_U16LE
    )
    imgData.create_dataset(
        "Fixed_Displacement_Values",
        data=values,
        dtype=H5T_IEEE_F32LE
    )

    # Close file and log completion
    file.close()
    logging.info('h5 export done!')

def mesh2voxelfe(mesh, templatefile, fileout, matprop=None, keywords=['NSET', 'ELSET'], eltype='C3D8', matpropbits=8, refnode=None, verbose=False):
    """Generate ABAQUS voxel Finite Element (FE) input file from 3D Unstructured Grid mesh data.
    The file written is an input file (.INP) in ABAQUS syntax that can be solved using ABAQUS or CALCULIX.
    The user can define a material mapping strategy for the conversion of local GVs to local material properties in the FE model.
    Material mapping laws are defined in separate template file(s) (see "prop.inp" and "property_temp_bone.inp" for examples).
    Boundary conditions, analysis type and output requests are defined in a separate template file (see "tmp.inp" for an example).
    Info on analysis definition at: https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-m-Sim-sb.htm#simacae-m-Sim-sb

    Parameters
    ----------
    mesh : meshio
        Unstructured grid mesh.
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
        * 'PROPERTY':
                Define an external material mapping law from template file. (Default = None)
                Use in combination with 'matprop' dictionary of material property files and corresponding GV ranges for the material mapping.
    eltype : str
        FE element type. The default is eight-node brick element (C3D8 and F3D8). See CalculiX node convention (Par. 6.2.1) at: http://www.dhondt.de/ccx_2.15.pdf
    matpropbits : int
        Bit depth for material mapping.
    refnode
        Reference node coordinates [REF_NODE_x, REF_NODE_y, REF_NODE_z] for kinematic coupling.
    verbose : bool
        Activate verbose output.
    """

    # verbose output
    if verbose:
        logging.basicConfig(level=logging.INFO)

    # variables initialization ##################################################
    # dictionaries (two) of GVs for each user-defined material property
    prop_GV = {}
    prop_GV2 = {}

    # dictionary of the GVs existing in the input model
    # existing_GV = []
    cell_GV = []
    GVmin = 0
    GVmax = 2 ** matpropbits - 1

    # check input model cell data
    try:
        mesh.cell_data
    except NameError:
        logging.warning('Mesh does not contain cell data. Material property mapping disabled.')
        binary_model = True
    else:
        binary_model = False
        keys = list(mesh.cell_data)
        cell_GV = mesh.cell_data.get(keys[0])
        if type(cell_GV) is list:
            cell_GV = cell_GV[0]
        # cell_GV = mesh.cell_data.get(keys[0])[0]
        logging.info('Found cell_data: {keyname}. cell_data range: {0} - {1}.'.format(min(cell_GV), max(cell_GV), keyname=keys[0]))
        if (min(cell_GV) == 0) & (max(cell_GV) == 1) & (issubclass(cell_GV.dtype.type, np.integer)):
            logging.warning('cell_data is binary. Material property mapping disabled.')
            binary_model = True

    # existing grey values
    existing_GV = np.unique(cell_GV)
    if max(existing_GV) > GVmax:
        GVmax = max(existing_GV)

    # initialize dictionaries of user material properties if given
    if binary_model:
        # GV range for binary dataset
        prop_GV[0] = (1, 1)
    else:
        if matprop is not None:
            # check matprop dictionary input
            if not "range" in matprop and len(matprop["file"]) > 1:
                raise IOError('Incorrect material property definition: GV ranges required when using multiple property files.')
            elif not "range" in matprop and len(matprop["file"]) == 1:
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
            logging.warning('matprop not given: material property mapping disabled.')

    # write ABAQUS *.inp output file ############################################
    logging.info('Start writing INP file')
    # open ABAQUS *.INP output file
    INP = open(fileout, 'w')

    # HEADER
    INP.write('** ---------------------------------------------------------\n')
    INP.write('** Abaqus .INP file written by Voxel_FEM script on {}\n'.format(datetime.now()))
    INP.write('** ---------------------------------------------------------\n')

    # NODE COORDINATES
    logging.info('Writing model nodes to INP file')
    INP.write('** Node coordinates from input model\n')
    INP.write('*NODE\n')

    # write only the existing nodes
    existing_nodes = np.unique(mesh.cells_dict["hexahedron"])
    for node in existing_nodes:
        INP.write('{0:10d}, {n[0]:12.6f}, {n[1]:12.6f}, {n[2]:12.6f}\n'.format(node, n=mesh.points[node]))

    # ELEMENTS AND ELEMENT SETS
    n_els = 0
    logging.info('Writing model elements to INP file')
    INP.write('** Elements and Element sets from input model\n')
    i = 0

    # for each existing model Grey Value
    for GV in existing_GV:
        # compose element SET definition string and write to output file
        elset_str = 'SET' + str(int(GV))
        INP.write('*ELEMENT, TYPE={0}, ELSET={1}\n'.format(eltype, elset_str))

        # for each cell with that Grey Value
        for cell in np.where(cell_GV == GV)[0]:
            # write cell index followed by a list of its nodes
            INP.write('{0},{n[0]},{n[1]},{n[2]},{n[3]},{n[4]},{n[5]},{n[6]},{n[7]}\n'.format(cell+1, n=mesh.cells[0][1][cell]))

    # NODE SETS
    if 'NSET' in keywords:
        INP.write('** Additional nset from voxel model. New generated nsets are:\n')
        INP.write('** {}, \n'.format(*list(mesh.point_sets)))

        # write node set string
        for nsetName in list(mesh.point_sets):
            INP.write('*NSET, NSET={0}\n'.format(nsetName))
            CR = 1 # carriage return
            # write node set indexes in lines of 10
            for node_i in mesh.point_sets[nsetName]:
                if CR < 10:
                    INP.write('{},'.format(node_i))
                else:
                    INP.write('{}\n'.format(node_i))
                    CR = 0
                CR = CR + 1
            INP.write('\n')
        logging.info('Additional nodes sets generated: {}'.format(list(mesh.point_sets)))

    # CELL SETS
    if 'ELSET' in keywords:
        INP.write('** Additional elset from voxel model. New generated elsets are:\n')
        INP.write('** {}, \n'.format(*list(mesh.cell_sets)))

        # write cell set string
        for csetName in list(mesh.cell_sets):
            INP.write('*ELSET, ELSET={}\n'.format(csetName))
            CR = 1 # carriage return

            # write cells set indexes in lines of 10
            for cell_i in mesh.cell_sets[csetName]:
                if CR < 10:
                    INP.write('{},'.format(cell_i))
                else:
                    INP.write('{}\n'.format(cell_i))
                    CR = 0
                CR = CR + 1
            INP.write('\n')
        logging.info('Additional cell sets generated: {}'.format(list(mesh.cell_sets)))

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
    logging.info('Model with {0} nodes and {1} elements written to file {fname}'.format(len(mesh.points), len(mesh.cells[0]), fname=fileout))

    return

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
