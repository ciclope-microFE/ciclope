from ciclope import recon_utils as ru
import numpy as np
from pyfabric import pyfabric

# I = ru.read_tiff_stack('/home/gianthk/Data/StefanFly_test/test_00__rec/test_00__rec_mini/recon_0000.tif')
# I = ru.read_tiff_stack('/home/gianthk/Data/2019.001.coop_TUberlin_simulierte_Mensch.iorig/trabecular_sample_mini2/2000L_crop_imgaussfilt_60micron_uint8_0000.tif')
# I = ru.read_tiff_stack('/home/gianthk/Data/StefanFly_test/test_00__rec/test_00__rec_mini2/recon_0000.tif')
I = ru.read_tiff_stack('/home/gianthk/Data/StefanFly_test/test_00__rec/test_00__rec/recon_00000.tiff')

pointset = np.array([[545, 542, 507], [649, 461, 522], [530, 301, 508]])

# fabric of one image
evecs, radii, evals, fabric_comp, DA = pyfabric.fabric(I[507 - 25:507 + 25, 542 - 25:542 + 25, 545 - 25:545 + 25])

print("here")