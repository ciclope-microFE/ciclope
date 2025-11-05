## Compute reaction forces from ParOSOl `.h5` model

```python
from ciclope.utils.postprocess import reaction_forces

output_file = '/home/giiori/terminus/external/tacosound/HR-pQCT_II/parosol/idx62.h5'

slice_level = 11
results = reaction_forces(output_file, slice_level)
```

### Results
---
- 'Z_min': 7.4538,
- 'Z_max': 7.7568,
- 'total_force': array([ -32.69463546,   32.95626426, -145.04631193]),
- 'F_tot': 152.29408112761442,
- 'num_nodes': 6078808,
- 'num_elements': 3950529,
- 'vs': 0.0303,
- 'nodes_z0_count': 154405,
- 'nodes_z1_count': 73097
---
