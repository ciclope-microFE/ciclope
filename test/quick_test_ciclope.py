import SimpleITK as sitk
import ciclope
import meshio
import ciclope.utils
from ciclope.utils.preprocess import remove_unconnected
import scipy
import numpy as np

import ciclope.utils.preprocess

# im = sitk.ReadImage("/Volumes/stamplab_terminus/processed/2025/compaCT/20250117_124458052387_crop_scaffold_masks/scaffold.mha")
im = sitk.ReadImage("/Users/gianthk/Desktop/scaffold.mha")

bw = remove_unconnected(sitk.GetArrayFromImage(im))
# ciclope.utils.recon_utils.plot_midplanes(bw)
bw2 = scipy.ndimage.binary_closing(bw,iterations=4, border_value=1)

# add plates
bw2 = ciclope.utils.preprocess.add_cap(np.transpose(bw2, axes=(1,2,0)), 10, 1)

# m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=[15.15e-6, 15.15e-6, 15.15e-6], meshtype='tetra', max_facet_distance=45e-6, max_cell_circumradius=45e-6)
m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=[0.0145, 0.0145, 0.0145], meshtype='tetra', max_facet_distance=0.06, max_cell_circumradius=0.06)

filename_mesh_out = '/Users/gianthk/Desktop/scaffold_tetramesh2.vtk'
m1.write(filename_mesh_out)

# m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=np.ones(3)*15.15e-6, max_facet_distance=1e-5, max_cell_circumradius=1e-5)
# m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=np.ones(3)*15.15e-6)

input_template = "/Users/gianthk/Code/ORMIR/ciclope/input_templates/tmp_comp_static_scaffold.inp"
filename_out = '/Users/gianthk/Desktop/scaffold_tetraFE.inp'
ciclope.core.tetraFE.mesh2tetrafe(m1, input_template, filename_out)