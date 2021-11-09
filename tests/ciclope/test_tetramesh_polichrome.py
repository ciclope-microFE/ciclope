from ciclope import ciclope
from ciclope import recon_utils as ru
from skimage import measure
from skimage.filters import threshold_otsu, gaussian
from scipy import ndimage

# resample factor
rf = 4

I = ru.read_tiff_stack('/home/gianthk/Data/TOMCAT/Kaya/D_single_h1h2_scale05/D_single_h1h2_scale050001.tif')
vs = [0.00325, 0.00325, 0.00325]


I = gaussian(I, sigma=1, preserve_range=True)

# resize the 3D data using spline interpolation of order 2
I = ndimage.zoom(I, 1/rf, output=None, order=2)

# correct voxelsize
vs = vs * rf

T = threshold_otsu(I)
BW = I > T
[labels, n_labels] = measure.label(~BW, None, True)
L = ciclope.remove_unconnected(BW)
labels[labels==1] = 0
L2 = ru.add_cap(L, cap_thickness=5, cap_val=2)
L3 = ru.add_cap(L2, cap_thickness=10, cap_val=3)

# viewer = napari.view_image(I)
# viewer.add_image(L3)

tetramesh = ciclope.cgal_mesh(L3, vs, 'tetra', 0.005, 0.01)
tetramesh.write('./test_tetra_polichrome.vtk')