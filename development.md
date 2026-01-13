## Development installation
Create and activate a virtual environment for development:
```commandline
conda env create -n ciclope
conda activate ciclope
```
Install [DXchange](https://dxchange.readthedocs.io/en/latest/index.html) (**Optional**, used in some of the examples) :
```shell
conda install -c conda-forge dxchange
```
Clone the git repository:
```shell
git clone https://github.com/ciclope-microFE/ciclope.git
```
Navigate to the repository folder and install the package using pip:
```shell
cd ciclope
pip install .
```
Or install the package locally for development with [Flit](https://flit.pypa.io/en/latest/index.html):
```commandline
flit install --symlink
```

### Basic installation (`voxelFE` only)
This section describes how to create a minimal environment for voxelFE only pipelines which don't require `CGAL`, `Eigen`, and `pygalmesh`.

```shell
micromamba create -n ciclope ipython python==3.11.11 jupyter numpy pymcubes scipy scikit-image meshio h5py -c conda-forge
git clone https://github.com/ciclope-microFE/ciclope.git
cd ciclope
pip3 install .
```

## Build instructions
Build package:
```commandline
python3 -m build
```
Upload to pypi:
```commandline
python3 -m twine upload --repository pypi dist/*
```
Uninstall the existing build:
```commandline
python -m pip uninstall ciclope
```

## HPC installation instructions
`ciclope` can be deployed on your High Performance Computing (HPC) environment. However, install instructions are environment specific and may vary a lot.
The following section shows how to install [`ParOSol`](https://github.com/reox/parosol-tu-wien) to run `ciclope` pipelines on a HPC cluster already containing modules for the required dependencies (`openmpi`, `eigen3`, `hdf5`, ...). 

Verify module availability and versions. The command `module spider` can be used to search for available software modules, showing all versions and their dependencies..
```bash
module avail openmpi
module key eigen
module key hdf5
module spider hdf5/1.8.23
```

Load the required modules:
```bash
module load openmpi/4.1.6
module load eigen/3.4.0
module load stack/2024-06 gcc/12.2.0
module load hdf5/1.8.23
module load cmake/3.27.7
```

Download and install [`ParOSol`](https://github.com/reox/parosol-tu-wien):
```bash
git clone https://github.com/reox/parosol-tu-wien.git
cd parosol-tu-wien/
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

Extract examples archive and run a test:
```bash
cd ..
tar xvf mesh.tar.bz2
mpirun -np 4 build/parosol mesh/h5/sphere.h5
```
