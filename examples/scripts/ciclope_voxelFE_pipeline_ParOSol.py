from skimage.filters import gaussian

from ciclope.utils.recon_utils import read_tiff_stack
from ciclope.utils.preprocess import remove_unconnected # , invert_images, convert_bmp_to_tiff, crop_and_resize_images
# from ciclope.utils.postprocess import calculate_total_force, circular_masks_BVTV, reaction_forces
from ciclope.core.voxelFE import vol2h5ParOSol

import numpy as np


# load data
input_file = './../../test_data/LHDL/3155_D_4_bc/cropped/3155_D_4_bc_0000.tif'
data_3D = read_tiff_stack(input_file)
vs = np.ones(3)*19.5e-3 # [mm]

# binarize
bw = remove_unconnected(gaussian(data_3D, sigma=1, preserve_range=True) > 63)

# ParOSol model generation
output_file = './../../test_data/LHDL/3155_D_4_bc/results/parosol/3155_D_4_bc.h5'

vol2h5ParOSol(
    bw, output_file,
    topDisplacement=-0.04, voxelsize=vs, poisson_ratio=0.3, young_modulus=18e3,
    topHorizontaFixedlDisplacement=True, locking_strategy="plane", plane_lock_num=3, verbose=False
)

# run ParOSol exe

# load reaction forces and calculate stiffness