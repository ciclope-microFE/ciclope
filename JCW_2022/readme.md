# biomechanics hackathon @ JCW_2022

## Anaconda and Jupyter crush course:

### Install conda env on your Jupyter
https://stackoverflow.com/questions/39604271/conda-environments-not-showing-up-in-jupyter-notebook
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
```

- Install packages from conda-forge channel
```commandline
conda install -c conda-forge dxchange
```

- Install packages with pip installer
```
pip install --upgrade PyMCubes
```