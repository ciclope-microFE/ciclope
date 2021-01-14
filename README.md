# CT2FE
#### (Computed Tomography to Finite Elements)

Utilities, scripts and notebooks for the conversion of 3D computed tomography data to finite element models. <br />
The folder [test_data](test_data) contains a test dataset for stress/displacement analysis. <br />
For further CalculiX examples visit [github.com/mkraska/CalculiX-Examples](https://github.com/mkraska/CalculiX-Examples) <br />

___

### Jupyter notebooks with CT2FE examples:
- [CT2FE_example01_voxelFE_static_Calculix.ipynb](CT2FE_example01_voxelFE_static_Calculix.ipynb)
![](test_data/example_01/masked_8bit_cap1.png)
    - [x] Load, inspect and write a 3D CT dataset
    - [x] Convert 3D stack to voxel-FE model for simulation in CalculX or Abaqus
        - [x] Local mapping of the dataset grey values to bone material properties
    - [x] Run simulation in Calculix
    - [x] Convert Calculix output to .VTK for visualization in Paraview
    - [x] Visualize simulation results in Paraview

- [CT2FE_example02_Slicer3Dmesher_Nlgeom_CalculiX.ipynb](CT2FE_example02_Slicer3Dmesher_Nlgeom_CalculiX.ipynb)
![](test_data/example_02/D_single_tens_Nlgeom.png)
    - [x] Load unstructured grid mesh preprocessed and generated with [3D Slicer](https://www.slicer.org/) <br />
    Tutorial for mesh generation in 3D Slicer using the -> [SlicerSegmentMesher module](https://github.com/lassoan/SlicerSegmentMesher#tutorial) 
    - [x] Read and modify mesh including:
        - [x] Material property definition
        - [x] Non-linear, quasi-static analysis definition: tensile test with material plasticity. Visit also [github.com/mkraska/CalculiX-Examples](https://github.com/mkraska/CalculiX-Examples/blob/master/Drahtbiegen/Zug/Zug.inp)
        - [x] Definition of boundary conditions
    - [x] Run simulation in Calculix
    - [x] Convert Calculix output to .VTK for visualization in Paraview
    - [x] Visualize simulation results in Paraview
___
### 2 DO:
**Pre-processing:**
- [x] add steel caps to the model
- [ ] 3D dataset embedding
- [ ] BMD-based material property calibration utilities

**Repo:**
- [ ] embedding
- [ ] BMD-based material property calibration utilities
- [ ] write midplanes images (.PNG)
- [ ] move TIFF (or other image formats) reading methods to separate library
- [ ] create module library

**Other:**
- [x] non-linear analysis (Nlgeom)
- [ ] non-linear analysis (bone damage [Rincón-Kohli 2008](10.1007/s10237-008-0128-z))
- [x] Calculix parallel (ccx_MT)
- [ ] BC and FE solution parameters definition as command line input (?)
- [ ] output different file format (.VTK)
- [ ] documentation pages (sphinx?)




