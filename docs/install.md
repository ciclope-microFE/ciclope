# Installation
For mesh generation, `ciclope` requires [pygalmesh-0.10.6](https://github.com/meshpro/pygalmesh), a Python frontend to [CGAL](https://www.cgal.org/).
Follow the [installation procedure](https://github.com/meshpro/pygalmesh#installation) for [CGAL](https://www.cgal.org/) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page).
After that, install pygalmesh with pip or conda:
```commandline
conda install -c conda-forge pygalmesh>=0.10.6
```
After installing pygalmesh, you can install `ciclope` using pip. The flag `[all]` will install optional dependencies needed to run full pipelines and examples.
```commandline
pip install ciclope[all]
```

Some examples will require [DXchange](https://dxchange.readthedocs.io/en/latest/index.html). You can install it with:
```shell
conda install -c conda-forge dxchange
```
---
## Testing
To verify your installation checkout this repository and run the tests with the command:
```commandline
cd test
python -m unittest -v test_ciclope.run_tests
```
---
## How to contribute
If you want to contribute to `ciclope` follow the installation steps below:

Create and activate a virtual environment for development:
```shell
conda env create -n ciclope
conda activate ciclope
```
Clone (or fork and clone) the git repository:
```shell
git clone https://github.com/gianthk/ciclope.git
```
Navigate to the repo folder and install the package using pip:
```shell
cd ciclope
pip install .
```
Or install the package locally for development with [Flit](https://flit.pypa.io/en/latest/index.html):
```shell
flit install --symlink
```
---
## Build instructions
Build package:
```shell
python -m build
```
Upload to pypi:
```shell
python -m twine upload --repository pypi dist/*
```
Uninstall the existing build:
```shell
python -m pip uninstall ciclope
```