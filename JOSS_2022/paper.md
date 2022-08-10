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

The python package `ciclope` processes micro Computed Tomography (microCT) data to generate Finite Element (FE) models.


# Statement of need
- microCT to finite elements pipelines are used to..
- our science case is musculoskeletal imaging but can be applied also to.. (examples)
- a variety of different pipelines was proposed `[@fernandez_nonlinear_2022; @megias_numerical_2022; @kox_heterogeneous_2022]`, none of which is fully open-source. As a result validation and comparison of results remains challenging.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from , and support from .. during the genesis of this project.

# References