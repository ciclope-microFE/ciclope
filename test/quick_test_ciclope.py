import SimpleITK as sitk
import ciclope
import meshio
from ciclope.utils.preprocess import remove_unconnected
import scipy
import numpy as np

# im = sitk.ReadImage("/Volumes/stamplab_terminus/processed/2025/compaCT/20250117_124458052387_crop_scaffold_masks/scaffold.mha")
im = sitk.ReadImage("/Users/gianthk/Desktop/scaffold.mha")

bw = remove_unconnected(sitk.GetArrayFromImage(im))
# ciclope.utils.recon_utils.plot_midplanes(bw)
bw2 = scipy.ndimage.binary_closing(bw,iterations=5, border_value=1)
# m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=[15.15e-6, 15.15e-6, 15.15e-6], meshtype='tetra', max_facet_distance=45e-6, max_cell_circumradius=45e-6)
m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=[15.15, 15.15, 15.15], meshtype='tetra', max_facet_distance=55.0, max_cell_circumradius=55.0)

filename_mesh_out = '/Users/gianthk/Desktop/scaffold_tetramesh.vtk'
# m1.write(filename_mesh_out)

# m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=np.ones(3)*15.15e-6, max_facet_distance=1e-5, max_cell_circumradius=1e-5)
# m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=np.ones(3)*15.15e-6)

input_template = "/Users/gianthk/Code/ORMIR/ciclope/input_templates/tmp_example01_comp_static_bone.inp"
filename_out = '/Users/gianthk/Desktop/scaffold_tetraFE.inp'
ciclope.core.tetraFE.mesh2tetrafe(m1, input_template, filename_out)