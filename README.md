# ciclope
Computed Tomography to Finite Elements.

[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/Naereen/StrapDown.js/blob/master/LICENSE)

[//]: # (more badges here: https://naereen.github.io/badges/)

**ciclope** processes micro Computed Tomography (microCT) data to generate Finite Element (FE) models. <br />

### Usage
**ciclope** can be run from the command line as a script. See the [examples](###Examples) below of this type of use. To view the help type
```commandline
python ciclope/ciclope.py -h
```

To use **ciclope** within python, import the module with
```python
from ciclope import ciclope
```

#### voxel-FE
![](test_data/trabecular_bone/trab_sample_mini3_UD3.png)
Read and segment a 3D dataset (TIFF stack) of trabecular bone. Generate **voxel-FE** model of linear elastic compression
```python
import numpy as np
from recon_utils import read_tiff_stack
from ciclope import ciclope

input_file = '/path/to/your/file.tiff'
input_template = "./../input_templates/tmp_example01_comp_static_bone.inp"

data_3D = read_tiff_stack(input_file)
vs = np.ones(3)*0.06 # [mm]

# segment and remove unconnected clusters
BW = data_3D > 142
L = ciclope.remove_unconnected(BW)

# generate CalculiX input file
ciclope.vol2voxelfe(L, input_template, 'test.inp', keywords=['NSET', 'ELSET'], voxelsize=vs)
```

#### tetrahedra-FE
![](test_data/steel_foam/B_matrix_tetraFE_mesh.png)
Read and segment 3D microCT dataset of steel foam sample
```python
import numpy as np
from skimage import measure, morphology
from recon_utils import read_tiff_stack
from ciclope import ciclope

input_file = '/your_path/steel_foam.tiff'

data_3D = read_tiff_stack(input_file)
vs = np.ones(3)*0.01625 # [mm]

# segment and remove unconnected clusters
BW = data_3D > 90
BW = morphology.closing(BW, morphology.cube(5))
L = ciclope.remove_unconnected(BW)
```

Generate mesh of tetrahedra with [pygalmesh](https://github.com/nschloe/pygalmesh)
```python
import pygalmesh
mesh = pygalmesh.generate_from_array(np.transpose(L, [2, 1, 0]).astype('uint8'), tuple(vs), max_facet_distance=0.02, max_cell_circumradius=0.1)
```
You can do the same within ciclope with
```python
mesh = ciclope.cgal_mesh(L, vs, 'tetra', 0.02, 0.1)
```

Generate **tetrahedra-FE** model of non-linear tensile test

```python
input_template = "./../input_templates/tmp_example02_tens_static_steel.inp"

# generate CalculiX input file
ciclope.mesh2tetrafe(mesh, input_template, 'test.inp', keywords=['NSET', 'ELSET'])
```

### ciclope pipeline 
The following table shows a general pipeline for FE model generation from CT data that can be executed with ciclope:

| # | Step | Description | **ciclope** flag |
|:-:|:-|:-|:-|
| 1. | **Load CT data** | | |
| 2. | **Pre-processing** | Gaussian smooth | `--smooth` |
| | | Resize image | `-r` |
| | | Add embedding | (not implemented yet) |
| | | Add caps | `--caps` |
| 3. | **Segmentation** | Uses Otsu method if left empty | `-t` |
| | | Remove unconnected voxels | |
| 4. | **Meshing** | Outer shell mesh of triangles | `--shell_mesh` |
| | | Volume mesh of tetrahedra | `--vol_mesh` |
| 5. | **FE model generation** | Apply Boundary Conditions | |
| | | Material mapping | `-m`, `--mapping` |
| | | Voxel FE | `--voxelfe` |
| | | Tetrahedra FE | `--tetrafe` |

---
### Notes on ciclope
* Tetrahedra meshes are generated with [pygalmesh](https://github.com/nschloe/pygalmesh) (a Python frontend to [CGAL](https://www.cgal.org/))
* High-resolution surface meshes for visualization are generated with the [PyMCubes](https://github.com/pmneila/PyMCubes) module.
* All mesh exports are performed with the [meshio](https://github.com/nschloe/meshio) module.
* **ciclope** handles the definition of material properties and FE analysis parameters (e.g. boundary conditions, simulation steps..) through separate template files. The folders [material_properties](/material_properties) and [input_templates](/input_templates) contain a library of template files that can be used to generate FE simulations.
  * Additional libraries of [CalculiX](https://github.com/calculix) examples and template files can be found [here](https://github.com/calculix/examples) and [here](https://github.com/calculix/mkraska)
___

### Examples
#### [Example #1 - voxel-FE model of trabecular bone](ciclope/ciclope_ex01_voxelFE_trabecularbone_CalculiX.ipynb) [![Made withJupyter](https://img.shields.io/badge/Made%20with-Jupyter-orange?style=for-the-badge&logo=Jupyter)](ciclope/ciclope_ex01_voxelFE_trabecularbone_CalculiX.ipynb)
![](test_data/trabecular_bone/U3.png)

The pipeline can be executed from the command line with:
```commandline
python ciclope/ciclope.py input.tif output.inp -vs 0.0606 0.0606 0.0606 --smooth -r 1.2 --vol_mesh --shell_mesh --voxelfe --template ./../input_templates/tmp_example01_comp_static_bone.inp -v
```

The example shows how to:
- [x] Load and inspect microCT volume data
- [x] Apply Gaussian smooth
- [x] Resample the dataset
- [x] Segment the bone tissue
- [x] Remove unconnected clusters of voxels
- [x] Convert the 3D binary to a voxel-FE model for simulation in CalculX or Abaqus
  - [x] Linear, static analysis; displacement-driven
  - [ ] Local material mapping (dataset grey values to bone material properties)
- [x] Launch simulation in Calculix
- [x] Convert Calculix output to .VTK for visualization in Paraview
- [x] Visualize simulation results in Paraview

#### [Example #2 - tetrahedra-FE model of stainless steel foam](ciclope/ciclope_ex02_tetraFE_steelfoam_CalculiX.ipynb) [![Made withJupyter](https://img.shields.io/badge/Made%20with-Jupyter-orange?style=for-the-badge&logo=Jupyter)](ciclope/ciclope_ex02_tetraFE_steelfoam_CalculiX.ipynb)
![](test_data/steel_foam/B_matrix_tetraFE_Nlgeom_results/PEEQ.gif)

The pipeline can be executed from the command line with:
```commandline
python ciclope/ciclope.py input.tif output.inp -vs 0.0065 0.0065 0.0065 --smooth -r 1.2 -t 90 --vol_mesh --tetrafe --template ./../input_templates/tmp_example02_tens_Nlgeom_steel.inp -v
```

The example shows how to:
- [x] Load and inspect synchrotron microCT volume data
- [x] Apply Gaussian smooth
- [x] Resample the dataset
- [x] Segment the steel
- [x] Remove unconnected clusters of voxels
- [x] Generate volume mesh of tetrahedra
- [x] Generate high-resolution mesh of triangles of the model outer shell (for visualization)
- [x] Convert the 3D binary to a tetrahedra-FE model for simulation in CalculX or Abaqus
  - [x] Non-linear, quasi-static analysis definition: tensile test with material plasticity. For more info visit: [github.com/mkraska/CalculiX-Examples](https://github.com/mkraska/CalculiX-Examples/blob/master/Drahtbiegen/Zug/Zug.inp)
  - [ ] Local material mapping
- [x] Launch simulation in Calculix
- [x] Convert Calculix output to .VTK for visualization in Paraview
- [x] Visualize simulation results in Paraview
 
___
### 2 DO:
- [ ] PACKAGING!

#### Documentation
- [X] link to [CalculiX examples](https://github.com/calculix/examples/tree/master/materials)
- [X] clean example notebooks
- [ ] add picture of CT 2 FE `ciclope.main()` pipeline
- [ ] API reference
- [ ] Usage
  - [X] as a script
  - [X] as a module - ciclope.methods
    - [X] `cgal_mesh`
    - [ ] `shell_mesh`
    - [X] `vol2voxelfe`
    - [X] `mesh2tetrafe`
  - [ ] prefect flow 

#### Pre-processing
- [x] add caps
- [X] write midplanes images (.PNG)
- [ ] 3D dataset embedding
- [ ] (line 598) custom copy of input_template with parameter substitution (driving node coordinates; displacement amplitude..)

#### pybonemorph
- [X] Center Of Mass
- [X] periosteum mask
- [ ] endosteum contour
- [ ] cortical bone mask

#### Post-processing
- [ ] CalculiX postprocessing
  - [X] process CalculiX output with [`dat2txt.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)
  - [X] produce Force displacement plot with `pandas`
  - [ ] cgx midplane plots
  - [ ] convergence plots with [`monitor.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)

#### Examples
- [X] CalculiX postprocessing
- [ ] steel caps
- [ ] embedding
- [ ] material property mapping with voxel-FE
- [ ] material property mapping with tetra-FE
  - [ ] multi-material tetrahedra-FE (femure with trabecular and cortical masks)
- [ ] **prefect** pipeline 1
  - [X] launch pipeline from master CSV table
  - [X] ciclope voxelFE pipeline
  - [ ] execute calculix
  - [ ] multiple load configurations
    - [ ] (comp, tens, shear, torsion, bending)
    - [ ] read input_template params from table (e.g. different displ amplitude)
  - [ ] postprocess results
    - [ ] midplanes Smises plots
    - [ ] generate report
    - [ ] write results to master



