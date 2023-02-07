# Usage
A general pipeline for FE model generation from microCT data is composed of the following steps
![FEM pipeline](./../examples/fem_pipeline.png)

To use **ciclope** within python, import the package with
```python
import ciclope
```
---
## Image pre-processing
`ciclope.utils` contains functions that help you read and pre-process 3D datasets for FE model generation.

Read 3D CT dataset stored as stack of TIFFs
```python
from ciclope.utils.recon_utils import read_tiff_stack

input_file = './test_data/LHDL/3155_D_4_bc/cropped/3155_D_4_bc_0000.tif'

data_3D = read_tiff_stack(input_file)
vs = np.ones(3) * 0.06  # voxelsize [mm]
```
`read_tiff_stack` reads all TIFF files (slices) contained in the `input_file` folder. The volume is stored in a 3D [`numpy.ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html) with size `[slices, rows, columns]`.

---
Segment the 3D image and remove unconnected voxels:
```python
from skimage import morphology
from ciclope.utils.preprocess import remove_unconnected

BW = data_3D > 142 # fixed global threshold
BW = morphology.closing(BW, morphology.ball(2)) # optional step
L = remove_unconnected(BW)
```
---
## Mesh and FE model generation
If you already have a mesh file, you can skip the mesh generation steps and use `ciclope` with 3D [`meshio`](https://github.com/nschloe/meshio) objects.
### voxel-FE
![](./../test_data/trabecular_bone/trab_sample_mini3_UD3.png)
Generate unstructured grid mesh of hexahedra (voxels)
```python
import ciclope
mesh = ciclope.core.voxelFE.vol2ugrid(data_3D, vs)
```

Generate CalculiX input file `.INP` for **voxel-FE** model of linear elastic compression test
```python
input_template = "./input_templates/tmp_example01_comp_static_bone.inp"
ciclope.core.voxelFE.mesh2voxelfe(mesh, input_template, 'foo.inp', keywords=['NSET', 'ELSET'])
```

### tetrahedra-FE
![](./../test_data/steel_foam/B_matrix_tetraFE_mesh.png)
Generate mesh of tetrahedra. `ciclope` uses [`pygalmesh`](https://github.com/nschloe/pygalmesh) for tetrahedra mesh generation
```python
mesh = ciclope.core.tetraFE.cgal_mesh(L, vs, 'tetra', max_facet_distance=0.2, max_cell_circumradius=0.1)
```

Generate CalculiX input file `.INP` for **tetrahedra-FE** model of non-linear tensile test
```python
input_template = "./input_templates/tmp_example02_tens_static_steel.inp"
ciclope.core.tetraFE.mesh2tetrafe(mesh, input_template, 'foo.inp', keywords=['NSET', 'ELSET'])
```
---
## Post-processing of results
`ciclope.utils.postprocess.paraviewplot` calls [ParaView](https://www.paraview.org/) to generate and save plots of a chosen model scalar field.

- Add path to your ParaView installation with
```python
import sys
sys.path.append('~/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib/python3.8/site-packages')
```

- Plot midplanes of the vertical displacement field `UD3`
```python
ciclope.utils.postprocess.paraview_plot('test_data/tooth/results/Tooth_3_scaled_2.vtk', slicenormal="xyz",
                                        RepresentationType="Surface", Crinkle=True, ColorBy=['U', 'D2'], Roll=90,
                                        ImageResolution=[1024, 1024], TransparentBackground=True,
                                        colormap='Cool to Warm')
```
| | | |
|:-------------------------:|:-------------------------:|:-------------------------:|
|![](./../test_data/tooth/results/Tooth_3_scaled_2_UD3_XY.png) | ![](./../test_data/tooth/results/Tooth_3_scaled_2_UD3_XZ.png) | ![](./../test_data/tooth/results/Tooth_3_scaled_2_UD3_YZ.png) |

- Plot midplanes of the Von Mises stress `S_Mises`
```python
ciclope.utils.postprocess.paraview_plot("test_data/tooth/results/Tooth_3_scaled_2.vtk", slicenormal="xyz",
                                        RepresentationType="Surface", Crinkle=False, ColorBy="S_Mises", Roll=90,
                                        ImageResolution=[1024, 1024])
```
| | | |
|:-------------------------:|:-------------------------:|:-------------------------:|
|![](./../test_data/tooth/results/Tooth_3_scaled_2_S_Mises_XY.png) | ![](./../test_data/tooth/results/Tooth_3_scaled_2_S_Mises_XZ.png) | ![](./../test_data/tooth/results/Tooth_3_scaled_2_S_Mises_YZ.png) |
See the Jupyter Notebooks in the [examples section](https://ciclope.readthedocs.io/en/latest/examples.html) for more examples of 3D data and results visualization.

---
## command-line use
**ciclope** pipelines can be run from the command line as a script. Scroll down and take a look at the [Examples](###Examples) for this type of use.
To view the command line script help run:
```shell
ciclope -h
```

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