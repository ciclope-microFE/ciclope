#### Development installation
Create and activate a virtual environment for development:
```commandline
conda env create -n ciclope
conda activate ciclope
```
Install [DXchange](https://dxchange.readthedocs.io/en/latest/index.html):
```shell
conda install -c conda-forge dxchange
```
Clone the git repository:
```shell
git clone https://github.com/gianthk/ciclope.git
```
Navigate to the repo folder and install the package using pip:
```shell
cd ciclope
pip install .
```
Or install the package locally for development with [Flit](https://flit.pypa.io/en/latest/index.html):
```commandline
flit install --symlink
```

#### Build instructions
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