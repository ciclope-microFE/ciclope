# pyfabric
Fabric tensor from 3D images.

pyfabric extracts the local fabric tensor from 3D image data.
pyfabric uses the 3D spatial Auto Correlation Function (ACF) to compute the local orientation and anisotropy of images.
The ACF is computed as:

$\int_\Omega \nabla u \cdot \nabla v~dx = \int_\Omega fv~dx$


To view the help type
```
python pyfabric.py -h
```

## to do:
- [X] fit ellipsoid to ACF
- [X] fabric_pointset
- [X] control ACF ROIsize
- [X] control ACF zoom
- [X] control ACF threshold
- [ ] `to01(ROIACF)>ACF_threshold` write with `numexpr`
- [ ] example prox femur
  - [X] peri mask
  - [X] peri mesh
  - [X] fabric pointset
  - [X] write fabric mesh
  - [X] Eigen-decomposition (get fabric tensor)
  - [ ] write Fabric tensor. Symmetric tensor components are expected to have the following order: XX, YY, ZZ, XY, YZ, XZ (see [here](https://kitware.github.io/paraview-docs/latest/python/paraview.simple.TensorGlyph.html))
  - [X] write anisotropy factor
- [ ] example demonstration ACF method
- [ ] test X, Y, Z orientations with known ellipsoid
- [ ] fabric_worm
- [ ] ROI checker (needed?)


