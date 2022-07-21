# ciclope 2 DO:
- [ ] voxel uFE with heterogeneous material mapping
- [ ] plot field mid-planes with Paraview
- [X] fix voxel and tetraFE models flippedlr
- [X] packaging
  - [X] allow usage of ciclope as a script and as module:
- [X] voxelFE: vol2ugrid (in meshio format!)
- [X] voxelFE: mesh2voxelFE

### Documentation
- [ ] comment FULL notebooks
- [ ] Read The Docs page
  - [ ] API reference
- [ ] add package to [jcmsk packages](https://jcmsk.github.io/packages.html)
- [ ] add E=f(GV) explanation and plotting
- [ ] add refs to material heterogeneity
- [ ] add picture of CT 2 FE `ciclope.main()` pipeline
- [ ] Usage
  - [X] as a script
  - [X] as a module - ciclope.methods
    - [X] `cgal_mesh`
    - [ ] `shell_mesh`
    - [X] `vol2ugrid` and `mesh2voxelfe`
    - [X] `mesh2tetrafe`
    - [ ] voxelFE: material mapping
    - [ ] tetraFE: material mapping
  - [ ] prefect flow 
- [X] Ex01 and Ex02 from JCW_2022 notebook pipelines (voxel and tetra)
- [X] add installation notes 2 readme
- [X] link to [CalculiX examples](https://github.com/calculix/examples/tree/master/materials)

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
- [ ] Field mid-planes plotting with PARAVIEW script
- [ ] CalculiX postprocessing
  - [X] process CalculiX output with [`dat2txt.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)
  - [X] produce Force displacement plot with `pandas`
  - [ ] cgx or [Paraview](https://www.paraview.org/Wiki/ParaView/Python/Screenshot) midplane plots
    ```shell
    export PYTHONPATH="/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib:/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib/python3.8/site-packages"```
  - [ ] convergence plots with [`monitor.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)

### Examples
- [ ] ex01: material property mapping with voxel-uFE
  - [ ] GV -> BMD -> E
- [ ] ex03: material property mapping with tetra-FE
  - [X] multi-material tetrahedra-FE (tooth+embedding)
  - [ ] plot displacement field midplanes with paraview
  - [ ] rescale U3 to reach 2 kN
  - [ ] (?) rotate tooth before embedding
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
- [X] CalculiX postprocessing
- [X] steel caps
- [X] embedding

### Validation
- [ ] check bone mass after material mapping vs ashing data
