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
# m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=np.ones(3)*15.15e-6, max_facet_distance=1e-5, max_cell_circumradius=1e-5)
m1 = ciclope.tetraFE.cgal_mesh(bw2, voxelsize=np.ones(3)*15.15e-6)
