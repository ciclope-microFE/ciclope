import SimpleITK as sitk
import ciclope
import meshio
import ciclope.utils
from ciclope.utils.preprocess import remove_unconnected
import scipy
import numpy as np
import mcubes
import pygalmesh
from ciclope.utils.recon_utils import read_tiff_stack
from skimage.filters import threshold_otsu, gaussian
from scipy import ndimage, misc

import ciclope.utils.preprocess

input_file = '/Users/gianthk/Code/ORMIR/ciclope/test_data/LHDL/3155_D_4_bc/cropped/3155_D_4_bc_0000.tif'
samplename = 'trab' # 'scaffold'
output_dir = '/Users/gianthk/Desktop/'

# read test image
data_3D = read_tiff_stack(input_file)
voxelsize = 0.0195 # [mm]
vs = np.ones(3)*voxelsize # [mm]

# gaussian filter
data_3D = gaussian(data_3D, sigma=1, preserve_range=True)

# downsample
resampling = 1

# resize the 3D data using spline interpolation of order 2
# data_3D = ndimage.zoom(data_3D, 1/resampling, output=None, order=2)

# correct voxelsize
vs = vs * resampling
voxelsize = voxelsize * resampling

# thresholding
bw = data_3D > 63 # from comparison with histology

# remove unconnected components
bw = remove_unconnected(bw)

# direct generation of tetrahedra mesh from array
filename_mesh_out = output_dir+samplename+'_tetramesh.vtk'
mesh_size_factor = 6
m1 = ciclope.tetraFE.cgal_mesh(bw, vs, 'tetra', max_facet_distance=min(vs), max_cell_circumradius=mesh_size_factor*min(vs))

# write mesh to file
m1.write(filename_mesh_out)

# generate mesh of hexahedra (voxel-mesh)
m3 = ciclope.core.voxelFE.vol2ugrid(bw, vs, verbose=True)

m3.cell_data['GV'] = m3.cell_data['GV'][0].astype('uint8')

# write mesh to file
filename_mesh_out = output_dir+samplename+'_voxelmesh.vtk'
m3.write(filename_mesh_out)

# pad binary image with a layer of zeros
pad_thickness = 3
bw = np.pad(bw, pad_thickness, 'constant', constant_values=0)

# generate surface mesh using marching cubes
# https://github.com/pmneila/PyMCubes/tree/master

vertices, triangles = mcubes.marching_cubes(bw, 0.5)

vertices = vertices * voxelsize
vertices = vertices[:,(2,1,0)] - pad_thickness*voxelsize # swap axes and remove padding

# write surface mesh
filename_surfacemesh_out = output_dir+samplename+'_surfacemesh.obj'
mcubes.export_obj(vertices, triangles, filename_surfacemesh_out)

# remeshing an existing surface mesh
# m11 = pygalmesh.remesh_surface(
#     filename_surfacemesh_out,
#     max_edge_size_at_feature_edges=4*voxelsize,
#     min_facet_angle=20,
#     max_radius_surface_delaunay_ball=7*voxelsize,
#     max_facet_distance=4*voxelsize,
#     verbose=True,
# )

# # write mesh to file
# filename_mesh_out = '/Users/gianthk/Desktop/'+samplename+'_surfacemesh_refined.vtk'
# m11.write(filename_mesh_out)

# volume mesh from surface mesh
m2 = pygalmesh.generate_volume_mesh_from_surface_mesh(
    filename_surfacemesh_out,
    min_facet_angle=20.0,
    max_radius_surface_delaunay_ball=mesh_size_factor*voxelsize,
    max_facet_distance=voxelsize,
    max_circumradius_edge_ratio=4.0,
    verbose=True,
    reorient=True,
)

# write mesh to file
filename_mesh_out = output_dir+samplename+'_tetramesh_from_surface.vtk'
m2.write(filename_mesh_out)
