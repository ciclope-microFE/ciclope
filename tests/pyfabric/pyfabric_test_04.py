import numpy as np
from scipy import ndimage, misc

from ciclope.recon_utils import read_tiff_stack, plot_midplanes
from pyfabric import pyfabric as pf
import meshio

input_file = '/home/gianthk/Data/2019.001.coop_TUberlin_simulierte_Mensch.iorig/2000L_crop_imgaussfilt_101micron_uint8/2000L_crop_imgaussfilt_101micron_uint8_0000.tif'

# Read 3D data
data_3D = read_tiff_stack(input_file)
vs = np.ones(3)*0.0101 # [mm]

filename_mesh_out = '/home/gianthk/PycharmProjects/CT2FE/test_data/trabecular_bone/2000L_crop_imgaussfilt_101micron_uint8.vtk'
mesh_trab = meshio.read(filename_mesh_out)

# cells barycenter coordinates
cells_bary = np.sum(mesh_trab.points[mesh_trab.cells[0][1][:]], 1)/mesh_trab.points[mesh_trab.cells[0][1][:]].shape[1]

evecs, radii = pf.fabric_pointset(data_3D, cells_bary/vs[0], ROIsize=50, ACF_threshold=0.33, ROIzoom=True, zoom_size=20, zoom_factor=2)



##########################
# plot the input image histogram
import matplotlib.pyplot as plt
fig2, ax2 = plt.subplots()
plt.hist(to01(ROIACF).ravel(), bins=100)
plt.show()