# Installation
Anaconda 2025.12 or newer is required.

1. Create a new conda environment
```commandline
conda create -n ciclope python=3.10
```

2. Activate the environment
```commandline
conda activate ciclope
```

3. Install IPyKernel and register the `ciclope` kernel to be seen by the main Jupyter installation
```commandline
conda install ipykernel
```

4. For mesh generation, `ciclope` requires [pygalmesh](https://github.com/meshpro/pygalmesh), a Python frontend to [CGAL](https://www.cgal.org/). Pygalmesh installed via (ana)conda from `conda-forge` channel already ships with a CGAL binary installation.
```commandline
conda install -c conda-forge pygalmesh
```

5. Install Ciclope
```
pip install ciclope[all]
```

[pygalmesh](https://github.com/meshpro/pygalmesh), a Python frontend to [CGAL](https://www.cgal.org/), should come with [CGAL](https://www.cgal.org/) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). If not, follow the [installation procedure](https://github.com/meshpro/pygalmesh#installation) for [CGAL](https://www.cgal.org/) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page).

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