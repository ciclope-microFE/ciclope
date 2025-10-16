import subprocess
import numpy as np

from skimage.filters import gaussian

from ciclope.utils.recon_utils import read_tiff_stack
from ciclope.utils.preprocess import remove_unconnected # , invert_images, convert_bmp_to_tiff, crop_and_resize_images
from ciclope.utils.postprocess import calculate_total_force, circular_masks_BVTV, reaction_forces
from ciclope.core.voxelFE import vol2h5ParOSol

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
# mpirun -np 8 build/parosol mesh/h5/sphere.h5

# Slurm call
# sbatch --ntasks=16 --time=2 --mem-per-cpu=100 --wrap="mpirun -np 16 build/parosol mesh/h5/sphere.h5"

# Define the sbatch command
command = [
    "sbatch",
    "--ntasks=16",
    "--time=2",
    "--mem-per-cpu=100",
    "--wrap=mpirun -np 16 build/parosol {output_file}"
]

# Run the command and capture the output
try:
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    print("Job submitted successfully!")
    print("SLURM output:", result.stdout)
except subprocess.CalledProcessError as e:
    print("Error submitting job:")
    print(e.stderr)

# see how to submit a joba rray here: https://www.gdc-docs.ethz.ch/EulerManual/site/jobs/

# load reaction forces and calculate stiffness
plane_lock_num = 10
slice_level = plane_lock_num + 1

results = reaction_forces(output_file, slice_level)

