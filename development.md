#### Development installation
Create and activate a virtual environment for development:
```commandline
conda env create -n ciclope
conda activate ciclope
```
Clone the git repository:
```commandline
git clone https://github.com/gianthk/ciclope.git
```
Navigate to the repo folder and install the package using pip:
```commandline
cd ciclope
pip install .
```
Or install the package locally for development with [Flit](https://flit.pypa.io/en/latest/index.html):
```commandline
flit install --symlink
```