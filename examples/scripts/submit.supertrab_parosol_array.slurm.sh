#!/bin/bash
#SBATCH --job-name=supertrab_parosol       #Name of the job
#SBATCH --array=1-80%10      #Array with 20 Jobs, always 10 running in parallel
#SBATCH --ntasks=1           #Requesting 1 node for each job (always 1)
#SBATCH --cpus-per-task=1    #Requesting 1 CPU for each job
#SBATCH --mem-per-cpu=2G     #Requesting 2 Gb memory per core and job
#SBATCH --time=4:00:00       #4 hours run-time per job
#SBATCH --output=supertrab_parosol_%a.log  #Log files


##########################################
echo "$(date) start ${SLURM_JOB_ID}"
##########################################

#Load the modules needed
# module load stack/2024-06 gcc/12.2.0
# module load hdf5/1.8.23 eigen cmake/3.27.7 python
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 eigen/3.4.0 hdf5/1.8.23

#define input and outputs
# OUT=SNPs
# if [ ! -e ${OUT} ]  ; then mkdir ${OUT} ; fi

##The internal variable of slurm (1-20 in our case; see header slurm) can be used to extract the names of the chromosomes. 
IDX=${SLURM_ARRAY_TASK_ID}
MODEL=$(sed -n ${IDX}p parosol_models_ds10_and_blur.lst)

#The ParOSol command
mpirun -np 1 ~/code/parosol-tu-wien/build/parosol ${MODEL}

##############################################
##Get a summary of the job
jeffrun -j ${SLURM_JOB_ID}
##############################################