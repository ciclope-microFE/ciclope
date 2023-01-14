# ciclope 2 DO:
- [X] plot field mid-planes with Paraview
- [ ] BCs as in [Cox et. al 2022. “Heterogeneous Tissue Modulus Improved Prediction of Mechanical Behavior in Osteoporotic Vertebral Cancellous Bone.”](https://www.biorxiv.org/content/10.1101/2021.11.30.470675v2)
- [X] voxel uFE with heterogeneous material mapping
- [X] fix voxel and tetraFE models flippedlr
- [X] packaging
  - [X] allow usage of ciclope as a script and as module:
- [X] voxelFE: vol2ugrid (in meshio format!)
- [X] voxelFE: mesh2voxelFE

### JOSS paper review
- [ ] rerun jupyter notebooks
- [ ] upload LHDL reproduction dataset to Zenodo or github
- [X] ethics LHDL dataset
- [X] single trabecula test dataset upload
**Tests:**
- [X] preprocess.remove_unconnected
- [X] voxelFE.vol2ugrid
- [X] voxelFE.mesh2voxelfe
- [X] tetraFE.shell_mesh
- [X] tetraFE.cgal_mesh
- [X] tetraFE.mesh2tetrafe
- [X] everything with `unittest`

### Documentation
- [X] comment FULL notebooks
- [X] Read The Docs page
  - [X] [`sphinx-build -nW --keep-going -b html docs/ docs/_build/html`](https://www.sphinx-doc.org/en/master/man/sphinx-build.html)
  - [X] import documentation
  - [ ] clean API reference
  - [X] clean notebooks Examples
  - [ ] usage
    - [ ] FE model generation - analysis definition
      - [ ] BCs
    - [ ] FE model generation - material mapping
      - [ ] add refs to material heterogeneity
      - [ ] add E=f(GV) explanation and plotting
    - [ ] `shell_mesh`
    - [X] as a module - ciclope.methods
      - [ ] voxelFE: material mapping
      - [ ] tetraFE: material mapping
      - [ ] post-processing of results
        - [ ] `ccx2paraview`
        - [X] PARAVIEW plot midplanes 
      - [X] `cgal_mesh`
      - [X] `vol2ugrid` and `mesh2voxelfe`
      - [X] `mesh2tetrafe`
    - [X] as a script
    - [ ] prefect flow 
- [X] add package to [ORMIR packages](https://ormircommunity.github.io/packages.html)
- [X] add picture of CT 2 FE `ciclope.main()` pipeline
- [X] Ex01 and Ex02 from JCW_2022 notebook pipelines (voxel and tetra)
- [X] add installation notes 2 readme
- [X] link to [CalculiX examples](https://github.com/calculix/examples/tree/master/materials)

### Pre-processing
- [x] add caps
- [X] write midplanes images (.PNG)
- [X] 3D dataset embedding
- [X] analysis template write with parameter substitution (driving node coordinates)
- [X] merge `pybonemorph` to `ciclope.utils.preprocess`
- [X] Center Of Mass
- [X] periosteum mask
- [ ] purge_mesh
- [ ] endosteum contour
- [ ] cortical bone mask

### Post-processing
- [X] Field mid-planes plotting with PARAVIEW
- [ ] CalculiX postprocessing
  - [X] process CalculiX output with [`dat2txt.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)
  - [X] produce Force displacement plot with `pandas`
  - [X] [Paraview](https://www.paraview.org/Wiki/ParaView/Python/Screenshot) midplane plots
    ```shell
    export PYTHONPATH="/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib:/home/gianthk/Applications/ParaView-5.9.0-RC1-MPI-Linux-Python3.8-64bit/lib/python3.8/site-packages"```
  - [ ] convergence plots with [`monitor.py`](https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts)

### Examples
- [x] ex01: material property mapping with voxel-uFE
  - [x] GV -> BMD -> E
- [X] ex03: material property mapping with tetra-FE
  - [X] plot displacement field midplanes with paraview
  - [X] rescale U3 to reach 2 kN
  - [X] comment NB
  - [X] multi-material tetrahedra-FE (tooth+embedding)
- [ ] Example **prefect** pipeline
  - [X] launch pipeline from master CSV table
  - [X] ciclope voxelFE pipeline
  - [ ] execute calculix
  - [ ] multiple load configurations
    - [ ] (comp, tens, shear, torsion, bending)
    - [ ] read input_template params from table (e.g. different displ amplitude)
- [ ] Example postprocessing of results
  - [ ] midplanes Smises plots with [Paraview](https://www.paraview.org/Wiki/ParaView/Python/Screenshot)
  - [ ] generate HTML report with [plotly](https://plotly.com/python/v3/html-reports/)
  - [ ] write results to master
  - [X] CalculiX postprocessing
- [X] steel caps
- [X] embedding

### Validation
- [ ] check bone mass after material mapping vs ashing data
