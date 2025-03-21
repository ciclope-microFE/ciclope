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
from ciclope.utils.recon_utils import bbox

def vol2ugrid(voldata, voxelsize=[1, 1, 1], GVmin=0, refnodes=False, verbose=False):
    """Generate unstructured grid mesh from 3D volume data.

    Parameters
    ----------
    voldata : ndarray
        3D voxel data.
    voxelsize : float
        3D model voxelsize.
    matpropbits : int
        Bit depth for material mapping.
    GVmin
        Minimum Grey Value considered for meshing. By default, GVmin=0: all zeroes are considered background.
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
    mesh : meshio
        Unstructured grid mesh.
    refnodes : list
        centroids on the model boundaries (X0, X1, Y0, Y1, Z0, Z1)
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
    # nodes arrays
    nodes = []
    existing_nodes = np.full((data_shape[0]+1)*(data_shape[1]+1)*(data_shape[2]+1), False)

    # cell nodes and cell grey value arrays
    cells = []
    cell_GV = []

    # dictionary of boundary node sets
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

    # reference nodes coordinates (one refnode per boundary)
    refnode_X0 = np.zeros(3)
    refnode_X1 = np.zeros(3)
    refnode_Y0 = np.zeros(3)
    refnode_Y1 = np.zeros(3)
    refnode_Z0 = np.zeros(3)
    refnode_Z1 = np.zeros(3)

    # dictionary of boundary cell sets
    cells_set = {
        'CELLS_X0': [],
        'CELLS_X1': [],
        'CELLS_Y0': [],
        'CELLS_Y1': [],
        'CELLS_Z0': [],
        'CELLS_Z1': []
    }

    # cell data ##################################################
    row_nodes = data_shape[2] + 1  # n nodes along y
    slice_nodes = (data_shape[2] + 1) * (data_shape[1] + 1)  # n nodes in 1 slice
    cell_i = 0

    logging.info('Calculating cell array')
    for slice in range(data_shape[0]):
        for row in range(data_shape[1]):
            for col in range(data_shape[2]):

                # get cell GV
                GV = voldata[(slice, row, col)]

                if GV > GVmin:
                    cell_i = cell_i + 1
                    # cell nodes indexes
                    # See eight-node brick cell (C3D8 and F3D8) node convention (Par. 6.2.1) at: http://www.dhondt.de/ccx_2.15.pdf
                    cell_nodes = [slice_nodes * slice + row_nodes * row + col,
                                  slice_nodes * slice + row_nodes * row + col + 1,
                                  slice_nodes * slice + row_nodes * (row + 1) + col + 1,
                                  slice_nodes * slice + row_nodes * (row + 1) + col,
                                  slice_nodes * (slice + 1) + row_nodes * row + col,
                                  slice_nodes * (slice + 1) + row_nodes * row + col + 1,
                                  slice_nodes * (slice + 1) + row_nodes * (row + 1) + col + 1,
                                  slice_nodes * (slice + 1) + row_nodes * (row + 1) + col]

                    # append to cell list
                    cells.append(cell_nodes)

                    # append to cell_GV list
                    cell_GV.append(GV)

                    # dictionary of cells belonging to boundary sets
                    if slice == 0:
                        cells_set['CELLS_Z0'].append(cell_i)
                    if slice == data_shape[0] - 1:
                        cells_set['CELLS_Z1'].append(cell_i)
                    if col == 0:
                        cells_set['CELLS_X0'].append(cell_i)
                    if col == data_shape[2] - 1:
                        cells_set['CELLS_X1'].append(cell_i)
                    if row == 0:
                        cells_set['CELLS_Y0'].append(cell_i)
                    if row == data_shape[1] - 1:
                        cells_set['CELLS_Y1'].append(cell_i)

                    # store if the cell nodes exist
                    for i in cell_nodes:
                        existing_nodes[i] = True

    # node coordinates and boundary node sets
    logging.info('Detecting node coordinates and boundary nodes')
    node_i = 0
    for slice in range(data_shape[0] + 1):
        for row in range(data_shape[1] + 1):
            for col in range(data_shape[2] + 1):
                # store node coordinates
                node_coors = [voxelsize[0] * col, voxelsize[1] * row, voxelsize[2] * slice]
                nodes.append(node_coors)

                if existing_nodes[node_i]:
                    # compose dictionary of boundary node sets
                    if slice == 0:
                        nodes_set['NODES_Z0'].append(node_i)
                        refnode_Z0 += node_coors
                    if slice == data_shape[0]:
                        nodes_set['NODES_Z1'].append(node_i)
                        refnode_Z1 += node_coors
                    if col == 0:
                        nodes_set['NODES_X0'].append(node_i)
                        refnode_X0 += node_coors
                    if col == data_shape[2]:
                        nodes_set['NODES_X1'].append(node_i)
                        refnode_X1 += node_coors
                    if row == 0:
                        nodes_set['NODES_Y0'].append(node_i)
                        refnode_Y0 += node_coors
                    if row == data_shape[1]:
                        nodes_set['NODES_Y1'].append(node_i)
                        refnode_Y1 += node_coors

                node_i = node_i + 1

    nodes_set['NODES_X0Y0Z0'] = nodes_set['NODES_Z0'][0:2]
    nodes_set['NODES_X0Y0Z1'] = nodes_set['NODES_Z1'][0:2]

    # reference nodes are the barycenters of boundary node sets
    refnode_X0 /= len(nodes_set['NODES_X0'])
    refnode_X1 /= len(nodes_set['NODES_X1'])
    refnode_Y0 /= len(nodes_set['NODES_Y0'])
    refnode_Y1 /= len(nodes_set['NODES_Y1'])
    refnode_Z0 /= len(nodes_set['NODES_Z0'])
    refnode_Z1 /= len(nodes_set['NODES_Z1'])

    refnodes_dict = {
                   "X0": refnode_X0,
                   "X1": refnode_X1,
                   "Y0": refnode_Y0,
                   "Y1": refnode_Y1,
                   "Z0": refnode_Z0,
                   "Z1": refnode_Z1,
    }

    # generate meshio object
    import meshio

    cells_dict = [("hexahedron", cells)]
    # mesh = meshio.Mesh(nodes, cells_dict, cell_data={"GV": [cell_GV]})
    mesh = meshio.Mesh(nodes, cells_dict, cell_data={"GV": [cell_GV]}, point_sets=nodes_set, cell_sets=cells_set)
    # mesh = meshio.Mesh(nodes, cells_dict, cell_data={"GV": [np.array(cell_GV, dtype='float')]})
    # mesh.write("foo.vtk")

    logging.info('Generated the following mesh with {0} nodes and {1} elements:'.format(len(mesh.points), len(mesh.cells[0])))
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
    
    #Calculating end using origin and dims 
    end=[origin[0]+dims[0], origin[1]+dims[1], origin[2]+dims[2]]
    
    #image dataset preparation
    logging.info('preparing Image dataset...')
    array_data = np.zeros((dims[0] + 1, dims[1] + 1, dims[2] + 1))
    
    lowBoundaryEls = []
    topBoundaryEls = []
    
    #logging.info('D shape ' + str(data_shape[0])+ ' ' + str(data_shape[1]) + ' ' + str(data_shape[2]))
    #logging.info('origin ' + str(origin[0])+ ' ' + str(origin[1]) + ' ' + str(origin[2]))
    #logging.info('dims ' + str(dims[0])+ ' ' + str(dims[1]) + ' ' + str(dims[2]))
        
    #fill output image data and get low and top boundary elements 
    # Define slices for the region of interest
    slice_start, slice_end = origin[0], end[0] + 1
    row_start, row_end = origin[1], end[1] + 1
    col_start, col_end = origin[2], end[2] + 1

    # Extract the sub-block from voldata
    sub_block = voldata[slice_start:slice_end, row_start:row_end, col_start:col_end]

    # Compute the new values in one operation
    result_block = sub_block * young_modulus

    # Assign directly to array_data
    array_data[:, :, :] = result_block
        
    # Find low boundary elements (slice == origin[0]) and non-zero values
    for i in range (0, plane_lock_num):
        low_boundary_mask = (sub_block[i, :, :] != 0)
        lowBoundaryEls2d = np.argwhere(low_boundary_mask)
        lowBoundaryEls2d = np.hstack((np.full((lowBoundaryEls2d.shape[0], 1), i), lowBoundaryEls2d))
        if (i == 0):
            lowBoundaryEls = lowBoundaryEls2d
        else:
            lowBoundaryEls = np.concatenate((lowBoundaryEls, lowBoundaryEls2d), axis=0)

    # Find top boundary elements (slice == end[0]) and non-zero values
    for i in range (0, plane_lock_num):
        top_boundary_mask = (sub_block[-(i+1), :, :] != 0)
        topBoundaryEls2d = np.argwhere(top_boundary_mask)
        topBoundaryEls2d = np.hstack((np.full((topBoundaryEls2d.shape[0], 1), slice_end - slice_start - (i+1)), topBoundaryEls2d))
        if (i == 0):
            topBoundaryEls = topBoundaryEls2d
        else:
            topBoundaryEls = np.concatenate((topBoundaryEls, topBoundaryEls2d), axis=0)
    
    #image dataset creation
    logging.info('creating Image dataset...')
    # Creazione di un array di esempio
    image=imgData.create_dataset("Image", (dims[0]+1,dims[1]+1,dims[2]+1), data=array_data, dtype=H5T_IEEE_F32LE)
    
    #setting voxel size
    logging.info('Setting voxel size')
    imgData.create_dataset("Voxelsize", data=voxelsize, dtype=H5T_IEEE_F64LE)

    #setting poisson ratio
    logging.info('Setting Poisson ratio')
    imgData.create_dataset("Poisons_ratio", data=poisson_ratio, dtype=H5T_IEEE_F64LE)
    
    #create nodes condition lists for top and low boundary els.
    lbNodesSet = set()
    tbNodesSet = set()

    #Creating boundary condition for nodes as sets 
    #each voxel has 8 nodes at the coordinates x and x+1, y and y+1, z and z+1
    #as reference read this 2d table:
        
    # Nodes:    0,0 --- 0,1 --- 0,2 --- 0,3 --- 0,4 --- 0,5 --- 0,6 
    # Elements:  |  0,0  |  0,1  |  0,2  |  0,3  |  0,4  |  0,5  |
    # Nodes:    0,0 --- 0,1 --- 0,2 --- 0,3 --- 0,4 --- 0,5 --- 0,6 
    # Elements:  |  0,0  |  0,1  |  0,2  |  0,3  |  0,4  |  0,5  |
    # Nodes:    0,0 --- 0,1 --- 0,2 --- 0,3 --- 0,4 --- 0,5 --- 0,6 
    # Elements:  |  0,0  |  0,1  |  0,2  |  0,3  |  0,4  |  0,5  |
    # Nodes:    0,0 --- 0,1 --- 0,2 --- 0,3 --- 0,4 --- 0,5 --- 0,6 
    
    #Voxels can have common nodes so we use a set to remove duplicates on nodes list 
    #the last column indicate the direction of the BCs
    
    for lbEl in lowBoundaryEls:
        lbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]  ,0) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]+1,0) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]  ,0) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]+1,0) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]  ,0) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]+1,0) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]  ,0) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]+1,0) )
        
        lbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]  ,1) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]+1,1) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]  ,1) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]+1,1) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]  ,1) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]+1,1) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]  ,1) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]+1,1) )

        lbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]  ,2) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]+1,2) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]  ,2) )
        lbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]+1,2) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]  ,2) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]+1,2) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]  ,2) )
        lbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]+1,2) )
    
    #if topHorizontaFixedlDisplacement is true X and Y displacements are fixed at the top
    if(topHorizontaFixedlDisplacement):    
        for tbEl in topBoundaryEls:
            tbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]  ,0) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]+1,0) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]  ,0) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]+1,0) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]  ,0) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]+1,0) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]  ,0) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]+1,0) )
            
            tbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]  ,1) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]+1,1) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]  ,1) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]+1,1) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]  ,1) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]+1,1) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]  ,1) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]+1,1) )

            tbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]  ,2) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]  , lbEl[2]+1,2) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]  ,2) )
            tbNodesSet.add( (lbEl[0]  , lbEl[1]+1, lbEl[2]+1,2) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]  ,2) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]  , lbEl[2]+1,2) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]  ,2) )
            tbNodesSet.add( (lbEl[0]+1, lbEl[1]+1, lbEl[2]+1,2) )
    else:
        for tbEl in topBoundaryEls:
            tbNodesSet.add( (tbEl[0]  , tbEl[1]  , tbEl[2]  ,2) )
            tbNodesSet.add( (tbEl[0]  , tbEl[1]  , tbEl[2]+1,2) )
            tbNodesSet.add( (tbEl[0]  , tbEl[1]+1, tbEl[2]  ,2) )
            tbNodesSet.add( (tbEl[0]  , tbEl[1]+1, tbEl[2]+1,2) )
            tbNodesSet.add( (tbEl[0]+1, tbEl[1]  , tbEl[2]  ,2) )
            tbNodesSet.add( (tbEl[0]+1, tbEl[1]  , tbEl[2]+1,2) )
            tbNodesSet.add( (tbEl[0]+1, tbEl[1]+1, tbEl[2]  ,2) )
            tbNodesSet.add( (tbEl[0]+1, tbEl[1]+1, tbEl[2]+1,2) )

	#Create Fixed_Displacement_Coordinates
    logging.info('Creating Fixed_Displacement_Coordinates')
    
    lbNsLen=len(lbNodesSet);
    tbNsLen=len(tbNodesSet);
    
    # Create a NumPy array of size lbNsLen+tbNsLen x 4
    fixDispCoord = np.zeros((lbNsLen+tbNsLen, 4), dtype=int)
    
    #put sets elements in the array
    lb_list = list(lbNodesSet)
    tb_list = list(tbNodesSet)
    fixDispCoord[:len(lb_list), :] = lb_list
    fixDispCoord[len(lb_list):, :] = tb_list
        
    fdc=imgData.create_dataset("Fixed_Displacement_Coordinates", data=fixDispCoord, dtype=H5T_STD_U16LE)
    
    #Create Fixed_Displacement_Values 
    logging.info('Creating Fixed_Displacement_Values')
    
    # Create a NumPy array of size lbNsLen+tbNsLen 
    fixDispValues = np.zeros((lbNsLen+tbNsLen), dtype=float)
            
    # All displacements are at 0 due to np.zeros so we need to add only top diplacement in the Z direction (ns[3] == 2)
    for idx, ns in enumerate(tbNodesSet):
        if(ns[3]==2):
            fixDispValues[lbNsLen+idx]=topDisplacement
        
    fdv=imgData.create_dataset("Fixed_Displacement_Values", data=fixDispValues, dtype=H5T_IEEE_F32LE)
    
    #finalizing 
    file.close()
    logging.info('hdf5 export done!')

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
