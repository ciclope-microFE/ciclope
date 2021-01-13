# CT2FE
(Computed Tomography to Finite Element)

This repository contains utilities, scripts and notebooks for the conversion of 3D computed tomography data to finite element models.
So far one example of a static stress/displacement analysis in Calculix or Abaqus is considered.

A dataset for running a test of stress/displacement analysis can be found [here]()

___

### Jupyter notebooks with CT2FE tests:
- [CT2FE_test01.ipynb](CT2FE_test01.ipynb)
    - Basic operations (load, inspect and write a 3D CT dataset)
    - Convert 3D stack to a voxel-FE model for simulation in Calculx or Abaqus
        - Mapping of the bone material properties based on local Grey Values is implemented
    - Run simulation in Calculix
    - Convert Calculix output to .VTK for visualization in Paraview
    - Visualize simulation results in Paraview

___
### 2 DO:
- Jupyter NOTEBOOK documenting one full example
    - [x] run complete test with Calculix
    - [ ] postprocessing Calculix output

- PREPROCESSING:
    - [x] add steel caps to the model
    - [ ] embedding
    - [ ] BMD-based material property calibration utilities

- FUNCTIONS:
    - [ ] embedding
    - [ ] BMD-based material property calibration utilities
    - [ ] write midplanes images (.PNG)
    - [ ] move TIFF (or other image formats) reading methods to separate library

- [ ] non-linear analysis (failure)
- [x] Calculix parallel (ccx_MT)
- [ ] BC and FE solution parameters definition as command line input (?)
- [ ] output different file format (.VTK)
- [ ] documentation pages (sphinx?)




