from src.ciclope import pybonemorph as pybm, recon_utils as ru
import pylab

# import cv2

# I = ru.read_tiff_stack('/home/gianthk/Data/StefanFly_test/test_00_/recon_phase/data_00100.tiff')
# I = tifffile.imread('/home/gianthk/Data/StefanFly_test/test_00_/recon_phase/data_00100.tiff')
I = ru.read_tiff_stack('/home/gianthk/Data/2019.001.coop_TUberlin_simulierte_Mensch.iorig/2000L_crop_imgaussfilt_101micron_uint8/2000L_crop_imgaussfilt_101micron_uint8_0000.tif')

# threshold
BW = I > 143

# test center of mass
# cx, cy = pybm.centerofmass(I[100,:,:]>31)
# cx, cy = com(BW)

# test perimask
peri = pybm.periosteummask(BW[640:680,:,:], 20)
# pylab.imshow(BW[655,:,:])
pylab.imshow(peri[20,:,:])
pylab.show()

print('here')