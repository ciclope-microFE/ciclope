# ciclope 2 DO:
- [X] voxelFE: vol2mesh (meshio standard)
- [X] voxelFE: mesh2voxelFE

### Documentation
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

### Pre-processing
- [x] add caps
- [X] write midplanes images (.PNG)
- [ ] purge_mesh
- [X] 3D dataset embedding
- [X] analysis template write with parameter substitution (driving node coordinates)

### pybonemorph
- [X] Center Of Mass
- [X] periosteum mask
- [ ] endosteum contour
- [ ] cortical bone mask

### Post-processing
- [ ] CalculiX postprocessing
  - [X] process CalculiX output with [`dat2txt.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)
  - [X] produce Force displacement plot with `pandas`
  - [ ] cgx or [Paraview](https://www.paraview.org/Wiki/ParaView/Python/Screenshot) midplane plots
    ```shell
    export PYTHONPATH="/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib:/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib/python3.8/site-packages"```
  - [ ] convergence plots with [`monitor.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)

### Examples
- [ ] ex03: material property mapping with tetra-FE
  - [ ] multi-material tetrahedra-FE (tooth+embedding)
  - [ ] rotate tooth before embedding
  - [ ] plot section of stress field with paraview
- [X] CalculiX postprocessing
- [X] steel caps
- [X] embedding
- [ ] material property mapping with voxel-FE
- [ ] **prefect** pipeline 1
  - [X] launch pipeline from master CSV table
  - [X] ciclope voxelFE pipeline
  - [ ] execute calculix
  - [ ] multiple load configurations
    - [ ] (comp, tens, shear, torsion, bending)
    - [ ] read input_template params from table (e.g. different displ amplitude)
  - [ ] postprocess results
    - [ ] midplanes Smises plots with [Paraview](https://www.paraview.org/Wiki/ParaView/Python/Screenshot)
    - [ ] generate HTML report with [plotly](https://plotly.com/python/v3/html-reports/)
    - [ ] write results to master
    