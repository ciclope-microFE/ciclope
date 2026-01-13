### Development installation
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

#### Basic installation (`voxelFE` only)
This section describes how to create a minimal environment for voxelFE only pipelines which don't require `CGAL`, `Eigen`, and `pygalmesh`.

```shell
micromamba create -n ciclope ipython python==3.11.11 jupyter numpy pymcubes scipy scikit-image meshio h5py -c conda-forge
git clone https://github.com/ciclope-microFE/ciclope.git
cd ciclope
pip3 install .
```
---
### Build instructions
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