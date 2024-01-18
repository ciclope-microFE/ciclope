#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ciclope module for beam Finite Element model generation
"""

import os
import logging
from datetime import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from scipy import ndimage, misc
from skimage.filters import threshold_otsu, gaussian
from skimage.morphology import skeletonize, skeletonize_3d,  medial_axis
from skimage.util import invert
from skimage import measure, morphology
from skimage import data
import sknw

import networkx as nx

from ciclope.utils.recon_utils import read_tiff_stack, plot_midplanes, plot_projections
from ciclope.utils.preprocess import remove_unconnected
import ciclope

def network_plot_3D(G, angle=30, save=False, point_size=20):
    """Generate 3D plot of graph data."""

    # Get node positions
    pos = nx.get_node_attributes(G, 'pts')

    # Get number of nodes
    n = G.number_of_nodes()

    # Get the maximum number of edges adjacent to a single node
    edge_max = max([G.degree(i) for i in range(n)])

    # Define color range proportional to number of edges adjacent to a single node
    colors = [plt.cm.plasma(G.degree(i) / edge_max) for i in range(n)]

    # 3D network plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Loop on the pos dictionary to extract the x,y,z coordinates of each node
    for key, value in pos.items():
        xi = value[0][0]
        yi = value[0][1]
        zi = value[0][2]

        # Scatter plot
        ax.scatter(xi, yi, zi, c=colors[key], s=point_size + point_size * G.degree(key), edgecolors='k', alpha=0.7)

    # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
    # Those two points are the extrema of the line to be plotted
    for i, j in enumerate(G.edges()):
        x = np.array((G.nodes[j[0]]['pts'][0][0], G.nodes[j[1]]['pts'][0][0]))
        y = np.array((G.nodes[j[0]]['pts'][0][1], G.nodes[j[1]]['pts'][0][1]))
        z = np.array((G.nodes[j[0]]['pts'][0][2], G.nodes[j[1]]['pts'][0][2]))
        # print(z)

        # Plot the connecting lines
        ax.plot(x, y, z, c='black', alpha=0.5)

    # Set the initial view
    ax.view_init(30, angle)

    # Hide the axes
    ax.set_axis_off()
    plt.show()

    return

def vol2graph(BWdata, vs=[1, 1, 1], verbose=False, thickness_method='max', ):

    # verbose output
    if verbose:
        logging.basicConfig(level=logging.INFO)

    # Distance transform from the image background
    logging.info('Computing distance transform')
    dist = ndimage.distance_transform_edt(BWdata)

    # 3D skeletonization with scikit.skeletonization_3d
    logging.info('Computing skeleton of given array')
    ske = skeletonize_3d(BWdata).astype(np.uint8)

    # Compute graph from the skeleton using the sknw package
    # build graph from skeleton
    logging.info('Calculating graph from skeleton')
    graph = sknw.build_sknw(ske)

    logging.info('Generated graph with {0} nodes and {1} elements:'.format(graph.number_of_nodes(), graph.number_of_elements()))

    # Assemble arrays for FE model creation
    #     nodes: [ID, x, y, z]
    #     elements: [ID_n1, ID_n2]
    #         NOTE: Node indexes start from 1!

    # initialize and assemble nodes array
    nodes = np.zeros((graph.number_of_nodes(),4), dtype=float)
    for n in range(0,graph.number_of_nodes()):
        nodes[n, 0]  = n+1
        nodes[n, 1:] = graph.nodes[n]['o']

    # elements array
    elements = np.array(list(graph.edges()))

    # initialize arc thickness array
    th = np.zeros((elements.shape[0]))

    logging.info('Calculating arc thickness')
    # get max thickness of the arc
    for n in range(0, len(th)):
        voxels = graph[elements[int(n),0]][elements[int(n),1]]['pts']
        th_array = np.zeros((voxels.shape[0],1))
        count = 0
        for v in voxels:
            th_array[count] = dist[v[0], v[1], v[2]]
            count += 1

        # assign trabecular thickness to the arc as twice the max distance from the arc skeleton points to the bone surface
        th[n] = 2*np.max(th_array)

    # Add mid-point of each element:
    logging.info('Adding arc mid-points')
    #     nodes: append mid-point node coordinates at the end of nodes array (n. of rows corresponiding to n. el will be added)
    #     elements: add third column corresponding to the row of the midpoint node in matrix nodes

    added_node_id = np.zeros((elements.shape[0],1))
    added_node_coors = np.zeros((elements.shape[0],4))

    # find mid-point node for each arc
    new_e_count = nodes.shape[0]
    for e in range(0,elements.shape[0]):
        n1 = elements[e,0]
        n2 = elements[e,1]

        coords = np.array([nodes[n1,1:],nodes[n2,1:]])
        n3 = np.mean(coords, axis=0)
        added_node_id[e] = new_e_count
        added_node_coors[e,0] = new_e_count+1
        added_node_coors[e,1:] = n3
        new_e_count += 1

    nodes_new = np.concatenate((nodes, added_node_coors), axis=0)
    elements_new = np.concatenate((elements, added_node_id), axis=1).astype('uint8')

    # Find Boundary nodes
    z_min = np.min(nodes_new[:,3])
    z_max = np.max(nodes_new[:,3])
    d_z = z_max-z_min
    clp=5/100
    print('Z min:', z_min)
    print('Z max:', z_max)

    idx_NODES_TOP = np.where(nodes_new[:,3] >= (z_min+(1-clp)*d_z))
    idx_NODES_BOTTOM = np.where(nodes_new[:,3] <= (z_min+clp*d_z))
    n_boundary = {'NODES_T': idx_NODES_TOP[0]+1,
                  'NODES_B': idx_NODES_BOTTOM[0]+1,
                  }

def graph2beamfe(nodes, elements, thickness, nsets_additional=None, fileout='pippo.inp', templatefile=None):
    """Generate ABAQUS beam Finite Element (FE) input file from graph data.
    The file written is an input file (.INP) in ABAQUS syntax that can be solved using ABAQUS or CALCULIX.
    Boundary conditions, analysis type and output requests are defined in a separate template file (see "tmp.inp" for an example).
    Info on analysis definition at: https://abaqus-docs.mit.edu/2017/English/SIMACAECAERefMap/simacae-m-Sim-sb.htm#simacae-m-Sim-sb

    Parameters
    ----------
    nodes : ndarry
        nodes array of type: [ID, x, y, z]
    elements : ndarray
        elements array of type: [ID_n1, ID_n2]
    thickness : ndarry
        element thickness array
    nsets_additional : dict
        Dictionary of boundary node sets as follows:
        nsets_additional = {
            'NODES_B': [0, 20, 32, ...],
            'NODES_T': [45, 52, 89, ...],
        }
    fileout : str
        Output .INP file.
    templatefile : str
        Analysis template file.
    eltype : str
        FE element type. The default is eight-node brick element (C3D8 and F3D8). See CalculiX node convention (Par. 6.2.1) at: http://www.dhondt.de/ccx_2.15.pdf
    verbose : bool
        Activate verbose output.
    """
    # open ABAQUS *.INP output file
    INP = open(fileout, 'w')

    # HEADER
    INP.write('** ---------------------------------------------------------\n')
    INP.write('** Abaqus .INP file written by ciclope')
    INP.write('** ---------------------------------------------------------\n')

    # NODE COORDINATES
    logging.info('Writing graph nodes to INP file')
    INP.write('** Node coordinates from input model\n')
    INP.write('*NODE, NSET=Nall\n')

    for node in nodes:
        INP.write('{0:4d}, {n[1]:12.6f}, {n[2]:12.6f}, {n[3]:12.6f}\n'.format(int(node[0]), n=node))

    # ELEMENTS AND ELEMENT SETS
    n_els = 0
    logging.info('Writing graph elements to INP file')
    INP.write('** Elements and Element sets from input model\n')
    INP.write('*ELEMENT, TYPE=b32, ELSET=SETALL\n')
    i = 0

    # for each element in the graph
    el_count = 1
    for el in elements:
        # write element index followed by a list of its nodes
        INP.write('{0:4d}, {n[0]:4d}, {n[2]:4d}, {n[1]:4d}\n'.format(el_count, n=el))
        el_count += 1

    # NODE SETS
    INP.write('** Additional nsets given:\n')
    if nsets_additional is not None:
        for key, arg in nsets_additional.items():
            INP.write('*NSET, NSET={0}\n'.format(key))
            for node in arg:
                INP.write('{}\n'.format(node))

            INP.write('\n')
            logging.info('Additional nset generated: {}'.format(key))

    # BOUNDARY CONDITIONS AND ANALYSIS DEFINITION
    # Open Abaqus analysis template file
    try:
        template = open(templatefile, 'r')
    except IOError('Abaqus template file {} not found.'.format(templatefile)):
        exit(1)
    logging.info('Reading Abaqus template file {}'.format(templatefile))

    # copy line by line info on model solution and boundary conditions from Abaqus analysis template file
    for line in template.readlines():
        # write line to output Abaqus file
        INP.write('{}'.format(line))

    template.close()
    INP.close()
    logging.info('Model with {0} nodes and {1} elements written to file {fname}'.format(len(nodes), len(elements), fname=fileout))

    return