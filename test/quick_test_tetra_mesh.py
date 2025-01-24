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
# input_file = '/Users/gianthk/Desktop/scaffold.mha'
# input_file = "/Volumes/stamplab_terminus/processed/2025/compaCT/20250117_124458052387_crop_scaffold_masks/scaffold.mha"
samplename = 'trab' # 'scaffold'

# read binary image
# data_3D = sitk.GetArrayFromImage(sitk.ReadImage(input_file))
# bw = data_3D

# read test image
data_3D = read_tiff_stack(input_file)
voxelsize = 19.5e-3
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
T = threshold_otsu(data_3D)
bw = data_3D > T

# remove unconnected components
bw = remove_unconnected(bw)

# simplify
# bw = scipy.ndimage.binary_closing(bw,iterations=1, border_value=1)
# bw = scipy.ndimage.binary_opening(bw,iterations=1, border_value=0)

# add endplates
# bw = ciclope.utils.preprocess.add_cap(np.transpose(bw, axes=(1,2,0)), 10, 1)

# generate mesh from array
filename_mesh_out = '/Users/gianthk/Desktop/'+samplename+'_tetramesh.vtk'
mesh_size_factor = 2
# m1 = ciclope.tetraFE.cgal_mesh(bw, vs, 'tetra', mesh_size_factor*min(vs), 2*mesh_size_factor*min(vs))

# write mesh to file
# m1.write(filename_mesh_out)

# generate microFE input file
# input_template = "/Users/gianthk/Code/ORMIR/ciclope/input_templates/tmp_comp_static_scaffold.inp"
# filename_out = '/Users/gianthk/Desktop/scaffold_tetraFE.inp'
# ciclope.core.tetraFE.mesh2tetrafe(m1, input_template, filename_out)

# pad with one voxel layer of zeros
bw = np.pad(bw, 3, 'constant', constant_values=0)

# smooth edges of the binary image
# bw = mcubes.smooth(bw)

# surface mesh using marching cubes
vertices, triangles = mcubes.marching_cubes(bw, 0.5)
vertices = vertices * voxelsize

# generate meshio mesh object
# cells = []

# for c in triangles:
#     cells.append(['triangle', c])

# m10 = meshio.Mesh(points=vertices, cells=cells)

# write surface mesh
filename_surfacemesh_out = '/Users/gianthk/Desktop/'+samplename+'_surfacemesh.obj'
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
    max_radius_surface_delaunay_ball=3*voxelsize,
    max_facet_distance=2*voxelsize,
    max_circumradius_edge_ratio=2*voxelsize,
    verbose=True,
    reorient=True,
)

# write mesh to file
filename_mesh_out = '/Users/gianthk/Desktop/'+samplename+'_tetramesh_from_surface.vtk'
m2.write(filename_mesh_out)
