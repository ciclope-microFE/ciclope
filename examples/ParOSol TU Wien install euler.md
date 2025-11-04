## bash command history 

```bash

module load stack/2024-05
module load gcc/13.2.0
module avail openmpi
module load openmpi/4.1.6
module key eigen3
module key eigen
module load eigen/3.4.0
module key hdf5
module load hdf5/1.14.5
module spider hdf5/1.14.5
module spider hdf5/1.12.2
module spider hdf5/1.8.23
module load hdf5/1.8.23
module spider hdf5/1.8.23
module load gcc/12.2.0
module load stack/2024-06 gcc/12.2.0
module load hdf5/1.8.23
git clone https://github.com/reox/parosol-tu-wien.git
cd parosol-tu-wien/
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release
module key cmake
module load cmake/3.27.7
cmake .. -DCMAKE_BUILD_TYPE=Release
make
cd ..
tar xvf mesh.tar.bz2
mpirun -np 46 build/parosol mesh/h5/sphere.h5
mpirun -np 16 build/parosol mesh/h5/sphere.h5
mpirun -np 8 build/parosol mesh/h5/sphere.h5
mpirun -np 4 build/parosol mesh/h5/sphere.h5

```

## launch jobs after installation
```bash
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 eigen/3.4.0 hdf5/1.8.23
sbatch --ntasks=16 --time=2 --mem-per-cpu=100 --wrap="mpirun -np 4 build/parosol mesh/h5/sphere.h5"
```
