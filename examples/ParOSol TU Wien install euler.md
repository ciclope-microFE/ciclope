```bash

module load stack/2024-06 gcc/12.2.0
module load hdf5/1.8.23 eigen cmake/3.27.7 python

sbatch --ntasks=16 --time=2 --mem-per-cpu=100 --wrap="mpirun -np 16 build/parosol mesh/h5/sphere.h5"
```
