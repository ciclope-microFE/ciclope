import recon_utils as ru
import pyfabric
import ellipsoid_fit as ef
import matplotlib.pyplot as plt

# I = ru.read_tiff_stack('/home/gianthk/Data/StefanFly_test/test_00__rec/test_00__rec_mini/recon_0000.tif')
# I = ru.read_tiff_stack('/home/gianthk/Data/2019.001.coop_TUberlin_simulierte_Mensch.iorig/trabecular_sample_mini2/2000L_crop_imgaussfilt_60micron_uint8_0000.tif')
I = ru.read_tiff_stack('/home/gianthk/Data/StefanFly_test/test_00__rec/test_00__rec_mini2/recon_0000.tif')

print(I.shape)

# get ACF
ACF = pyfabric.ACF(I)

# zoom center
ACF = pyfabric.zoom_center(ACF)

# convert to 0-1 range
ACF = ru.to01(ACF)

# binarize through FWHM
bw = ACF > 0.5

# get coordinates of envelope
coors = pyfabric.envelope(bw)

# plot envelope points
ax = pyfabric.scatter_plot(coors)

# fit ellipsoid
center, evecs, radii, v = ef.ellipsoid_fit(coors)

# plot
ef.ellipsoid_plot(center, radii, evecs, ax=ax, plot_axes=True, cage_color='orange')
plt.show()

# # plotting
# # # ru.plot_midplanes(ACF)
# pylab.imshow(ACF[int(I.shape[0]/2),:,:])
# pylab.show()
#
# pylab.imshow(ACF[:,int(I.shape[1]/2),:])
# pylab.show()
#
# pylab.imshow(ACF[:,:,int(I.shape[2]/2)])
# pylab.show()