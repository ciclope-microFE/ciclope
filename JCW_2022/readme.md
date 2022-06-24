# Biomechanics group hackathon @ 2022_JCW
#### Building the Jupyter Community in MSK Imaging Research - A Jupyter Community Workshop sponsored by [NUMFocus](https://numfocus.org/)
#### Organized by and for the [Jupyter Community in MSK Imaging Research](https://jcmsk.github.io/) 
Dates: June 9-11, 2022

---
### Group activities presentation
- Find [here](QMSKI2022-ORMRWorkshop-Biomechanics.pdf) a presentation of the activities and report from the biomechanics group @ 2022_JCW 
### Hackathon Jupyter notebooks
- [CT to voxel uFE pipeline](notebooks/JCW_voxeluFE_CalculiX.ipynb)
- [CT to tetra uFE pipeline](notebooks/JCW_tetrauFE_CalculiX.ipynb)

---
## Anaconda and Jupyter crush course:

### Install conda env on your Jupyter
Source [here](https://stackoverflow.com/questions/39604271/conda-environments-not-showing-up-in-jupyter-notebook)
```commandline
source activate ciclope
python -m ipykernel install --user --name ciclope --display-name "conda (ciclope)"
```

- Link to conda executable from within console
```commandline
where conda
path %PATH%;C:\ProgramData\Anaconda3\Scripts\
```

- Install packages
``` 
conda install numpy
conda install matplotlib
conda install scikit-image
conda install tqdm
```

- Install packages from conda-forge channel
```commandline
conda install -c conda-forge dxchange
conda install -c conda-forge meshio=5.0.0
conda install -c conda-forge pygalmesh
conda install -c conda-forge ccx2paraview
```

- Install packages with pip installer
```
pip install --upgrade PyMCubes
pip install --upgrade pypng
pip intsall itkwidgets
pip install vtk
```

