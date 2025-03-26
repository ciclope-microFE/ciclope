import os
import numpy as np

from ciclope.utils.recon_utils import read_tiff_stack, plot_midplanes
import ciclope.utils.preprocess
from ciclope.utils.preprocess import remove_unconnected
import ciclope

data_folder = '/Volumes/stamplab_terminus/processed/Zeiss730/2025/02/01_Scaffold-nocell_2025-02-07_111220/masks/'
datasets = ['scaffold2_clean', 'scaffold3_clean'] # 'scaffold1_clean', 
output_folder = '/Volumes/stamplab_terminus/processed/Zeiss730/2025/02/01_Scaffold-nocell_2025-02-07_111220/meshes/' 

for dataset in datasets:

    input_file = os.path.join(data_folder, dataset, 'slice_001.tiff')
    data_3D = read_tiff_stack(input_file)
    vs = np.ones(3)*9.4e-3 # [mm]

    # invert the image
    BW = np.invert(data_3D)

    # crop the image along x and y
    crop_size_z = 10
    crop_size_y = 12
    BW2 = BW[crop_size_z:-crop_size_z, crop_size_y:-crop_size_y, :]

    # remove unconnected components
    BW3 = remove_unconnected(BW2)

    # pad image along X
    pad_width = 40

    # Create a padding array with True values
    padding_shape = (BW3.shape[0], BW3.shape[1], pad_width)
    BW4 = np.ones(padding_shape, dtype=bool)

    # Concatenate the original array with the padding array along the 3rd dimension
    BW4 = np.concatenate((BW4, BW3.astype('bool'), BW4), axis=2, dtype=bool)

    # create mesh of tetrhedra
    max_facet_distance=0.02
    max_cell_circumradius=0.05

    m1 = ciclope.tetraFE.cgal_mesh(BW4, voxelsize=vs, meshtype='tetra', max_facet_distance=max_facet_distance, max_cell_circumradius=max_cell_circumradius)

    # write mesh output
    filename_mesh_out = os.path.join(output_folder, dataset+'_VoidMesh_MaxFacetDistance'+str(max_facet_distance)+'_MaxCellCircradius'+str(max_cell_circumradius)+'.vtk')
    m1.write(filename_mesh_out)
    filename_mesh_out = os.path.join(output_folder, dataset+'_VoidMesh_MaxFacetDistance'+str(max_facet_distance)+'_MaxCellCircradius'+str(max_cell_circumradius)+'.nas')
    m1.write(filename_mesh_out)
    # filename_mesh_out = os.path.join(output_folder, datasets[0]+'voidmesh.stl')
    # m1.write(filename_mesh_out)

