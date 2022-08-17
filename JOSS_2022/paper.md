---
title: 'Ciclope: micro Computed Tomography to Finite Elements'

tags:
  - Python
  - computed tomography
  - finite elements
  - image processing
  - engineering
  - simulation
  - biomechanics

authors:
  - name: Gianluca Iori
    orcid: 0000-0001-8611-3281
    affiliation: 1
    corresponding: true
  - name: Fulvia Taddei
    orcid: 0000-0002-8342-7434
    affiliation: 2
  - name: Enrico Schileo
    orcid: 0000-0003-3515-1231
    affiliation: 2
  - name: Gianluigi Crimi
    affiliation: 2
  - name: Giulia Fraterrigo
    affiliation: 2
  - name: Martino Pani
    orcid: 0000-0002-5786-4462
    affiliation: 3

affiliations:
  - name: Synchrotron-light for Experimental Science and Applications in the Middle East, Jordan
    index: 1
  - name: Rizzoli Orthopaedic Institute, Bologna, Italy
    index: 2
  - name: School of Mechanical and Design Engineering, University of Portsmouth, UK
    index: 3

date: 15 August 2022
bibliography: pippo.bib

---

# Summary
We present the python package `ciclope`, which processes micro Computed Tomography (microCT) data to generate Finite Element (FE) models. 

# Statement of need
Tissue level micro-Finite Element (microFE) models derived from laboratory or synchrotron microCT images can provide non-destructive assessments of the bone mechanical properties. The technique is used to investigate the effect of pathologies, treatment and remodelling on the mechanical response of bone at the tissue level, and is applied both to human and animal samples. Linear elastic microFE simulations are implemented to back-calculate the tissue elastic modulus (Bayraktar et al. 2004), understand deformation mechanisms (Zauel et al. 2005), or even predict failure (Pistoia et al. 2002) of cancellous bone, as well as to estimate the stiffness of whole bones from small animals (oliviero). Several groups highlighted the importance of an accurate description of boundary conditions and of validating model predictions with experimental measurements (e.g. with Digital Volume Correlation). Nevertheless, the verification of microFE procedures is hindered by: (i) a general absence of standardized, open procedures, (ii) the number of different, non-open source FE solvers and (iii) limited access to open datasets and to X-ray imaging and high-performance computing facilities.

Different pipelines for the generation of microFE models of trabecular bone have been proposed [@fernandez_nonlinear_2022; @megias_numerical_2022; @cox_heterogeneous_2022], none of which is fully open-source. As a result validation and comparison of results remains challenging.

The development of open-source and reproducible microFE workflows is expected to facilitate and support the validation of biomechanical studies, strengthening at the same time the synergy with other fields of microFE application such as concrete, fiber composites and porous materials research.

We present a fully-open-source pipeline from microCT data preprocessing to FE model generation, solution and postprocessing (visualization) of results.

Our science case is musculoskeletal imaging, but can microCT-derived FE models be applied to various fields (add examples)

# Design
**Ciclope** is composed of a core library of modules for FE model generation (`ciclope.core`), and a library of utilities for image and FE model pre- and post-processing (==ciclope.utils==) that can be imported and used within Python. Additionally, the `ciclope.py` script generated during package installation <mark>`ciclope.py`</mark> allows to launch microCT-to-FE pipelines directly from the commandline.

![Design of ciclope, and application to a pipeline for FE model generation from microCT data.\label{fig:design}](./../docs/ciclope_design.png){width=100%}

A pipeline for the generation and solution of a FE model derived from 3D microCT data is shown in the central part of figure \autoref{fig:design}.
**Image pre-processing**: a microCT dataset is loaded in python and segmented to isolate bone voxels and background. A connectivity check is performed to remove spurious, unconnected structures. The 3D image can be smoothed, rotated, cropped and resampled to lower resolution. Embedding layers and steel caps can be added to simulate experimental conditions of mechanical testing.
**Meshing**: **ciclope** allows to create several types of FE meshes. Image voxels can be directly converted to 8-node, hexahedral brick elements with the `voxelFE.py` module. Alternatively, meshes of 4-node tetrahedra can be generated with [@cgal] (`tetraFE.py` module). Finally, the `trussFE.py` module allows to generate a mesh of 2-node beam elements, where each beam represents a single trabecula, and has a local trabecular thickness associated to it.
**FE model generation**: the mesh is converted to an `.INP` input file for the FE solver. Within this process, the user can define the model material properties and the type of FE analysis (i.e. boundary conditions, analysis type and steps, requested outputs) through separate `.TMP` template files. Libraries of `material_properties` and `input_templates` are provided. Additional resources can be found online [here](https://github.com/calculix/examples) and [here](https://github.com/calculix/mkraska). For voxel-FE model generation, different **material mapping** strategies can be used: uniform tissue material properties (elastic modulus and poisson ratio) can be applied to all bone voxels. Alternatively, the voxel grey values (GV) can be converted to heterogeneous material properties using a mapping law defined by the user.\

[comment]: <> (bone mineral density BMD through a calibration rule obtained scanning a hydroxyapatite phantom. After this, an empirical law is used to convert local BMD to tissue elastic moduli Bourne_2004; garcia_2008.)
[comment]: <> (The pipeline is composed of the following steps:)
[comment]: <> (1. **microCT image pre-processing**: after reading in python a microCT dataset, the 3D volume can be cropped and aligned according to the desired direction of load, smoothed to remove noise with a Gaussian kernel, and resampled to lower image resolution. A binary mask of the bone tissue is generated thresholding bone voxels. Several global Otsu; Ridler_1978, or local adaptive thresholding ,..., techniques have been proposed Kim 2006. Embedding layers and steel caps can be added to simulate the experimental conditions of mechanical testing.) 

## The ciclope ecosystem
**Ciclope** relies on several dependencies:

- Voxel and tetrahedra mesh exports are performed with [meshio](https://github.com/nschloe/meshio) [@meshio].
- Tetrahedra meshes are generated with [pygalmesh](https://github.com/nschloe/pygalmesh), a Python frontend to [CGAL](https://www.cgal.org/) [@pygalmesh].
- High-resolution surface meshes for visualization are generated with [PyMCubes](https://github.com/pmneila/PyMCubes).
- The FE `.INP` files generated by **ciclope** can be solved using the free software [CalculiX](https://github.com/calculix) [@calculix] or [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/).
- For the visualization of FE results, ciclope relies on [ParaView](https://www.paraview.org/) [@paraview] and the CalculiX to ParaView converter [`ccx2paraview`](https://github.com/calculix/ccx2paraview) [@ccx2paraview].

[comment]: <> (Dxchange @decarlo_2014)

# Examples
`ciclope` was initially created to solve the problem of brittle crack propagation using thermodynamically consistent framework [@kaczmarczyk2017energy]. 
Over time, the domain of applications expanded to include computational homogenisation [@ullah2019unified], bone remodelling and fracture [@lew2020numerical], modelling of gel rheology [@richardson2018multiphysics] and acoustics problems. Moreover, `MoFEM` includes an extensive library of example applications such as soap film, solid shell, topology optimisation, phase field fracture, Navier-Stokes flow, cell traction microscopy, bone remodelling, configurational fracture, plasticity, mortar contact, magnetostatics and acoustic wave propagation as shown in \autoref{fig:examples}.

![Examples of microCT data to microFE model generation with `ciclope`.\label{fig:examples}](./../docs/nosignal.jpg){width=100%}

# Conclusions
If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements
We acknowledge contributions from , and support from .. during the genesis of this project.
The **ciclope** project was partially developed during the Jupyter Community Workshop [“Building the Jupyter Community in Musculoskeletal Imaging Research”](https://github.com/JCMSK/2022_JCW) sponsored by [NUMFocus](https://numfocus.org/).

# References