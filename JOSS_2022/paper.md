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
  - name: Martino Pani
    orcid: 0000-0002-5786-4462
    affiliation: 2
affiliations:
  - name: Synchrotron-light for Experimental Science and Applications in the Middle East, Jordan
    index: 1
  - name: School of Mechanical and Design Engineering, University of Portsmouth, UK
    index: 2
date: 05 August 2022
bibliography: pippo.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
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
![Design of ciclope, and application to a common pipeline for FE model generation from microCT data. \label{fig:design}](./../docs/ciclope_design.png)
The central diagram of figure \autoref{fig:design} illustrates a pipeline for the generation, solution, and postprocessing of results from a microFE model derived from 3D, microCT data.
The pipeline is composed of the following steps:
1. **microCT image pre-processing**: after reading in python a microCT dataset, the 3D volume can be cropped and aligned according to the desired direction of load, smoothed to remove noise with a Gaussian kernel, and resampled to lower image resolution. A binary mask of the bone tissue is generated thresholding bone voxels. Several global (Otsu; Ridler_1978) or local adaptive thresholding (refs) techniques have been proposed (Kim 2006). Embedding layers and steel caps can be added to simulate the experimental conditions of mechanical testing. 
2. **mesh generation**:
3. **FE model generation**:
   3.1. **Material mapping**: Uniform tissue material properties (elastic modulus and poisson ratio) are generally applied to all bone voxels. Alternatively, voxel grey values (GV) can be converted to bone mineral density (BMD) through a calibration rule obtained scanning a hydroxyapatite phantom. After this, an empirical law is used to convert local BMD to tissue elastic moduli (Bourne_2004; garcia_2008). 
4. **FE model solution**:

(ii) a core module for the generation of a microFE model from the pre-processed dataset, (iii) model solution with an open-source FE solver and (iv) post-processing of results and 3D visualization routines. 
The python script stack2Abaqus.py converts the preprocessed 3D dataset to a microFE input model file for solution with the open-source code CalculiX or with the commercial software ABAQUS. Image voxels are directly converted to 8-node hexahedral brick elements. The user can define either constant material properties or separate material mapping laws for different GV ranges. In this way, different material properties can be applied, for example, to bone tissue and steel or embedding material. Automatically generated boundary node sets are used to prescribe a uniaxial compression test up to 1% strain. The model is solved with a multi threaded version of CalculiX. The sample stiffness is calculated from the displacements and reaction forces written by CalculiX in separate output files. The visualization of the deformed 3D structure in Paraview is also illustrated after conversion with ccx2paraview.

## Ciclope modules
The package is composed of a core module containing methods for FE model generation (`ciclope.core`), and a module of utilities for image and FE model pre- and post-processing (`ciclope.utils`).

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

## The ciclope ecosystem
**Ciclope** requires several dependencies. The following is a list of the main external packages required for FE model generation and solution:
* All mesh exports (voxel and tetrahedra Finite Elements) are performed with the [meshio](https://github.com/nschloe/meshio) module.
* Tetrahedra meshes are generated with [pygalmesh](https://github.com/nschloe/pygalmesh) (a Python frontend to [CGAL](https://www.cgal.org/)).
* High-resolution surface meshes for visualization are generated with the [PyMCubes](https://github.com/pmneila/PyMCubes) module.
* **Ciclope** generates Finite Element `.INP` files that can be solved using both [CalculiX](https://github.com/calculix) and [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/).
* The definition of material properties and of the FE analysis parameters (e.g. boundary conditions, simulation steps..) is handled through separate template files. The folders [material_properties](https://github.com/gianthk/ciclope/tree/master/material_properties) and [input_templates](https://github.com/gianthk/ciclope/tree/master/input_templates) contain a library of template files that can be used to generate FE simulations.
  * Additional libraries of [CalculiX](https://github.com/calculix) examples and template files can be found [here](https://github.com/calculix/examples) and [here](https://github.com/calculix/mkraska)
* For the post-processing of FE results, **ciclope** uses [ParaView](https://www.paraview.org/) and the CalculiX to ParaView converter [`ccx2paraview`](https://github.com/calculix/ccx2paraview).

## Usage
Both modules can be imported and used within Python.
Pipelines of FE model generation can be launched from the commandline using the `ciclope.py` script generated during installation. See the section [usage](usage) for more details.

# Examples
description of examples..

# Conclusions
If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

We acknowledge contributions from , and support from .. during the genesis of this project.

# References