#!/bin/bash -x
#SBATCH --account=rcongm
#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --output=./slurm_dales-%j.out
#SBATCH --error=./slurm_error-%j.out

# *** start of job script ***
# Note: The current working directory at this point is
# the directory where sbatch was executed.

 module load Stages/2024
 module load GCC/12.3.0   ParaStationMPI/5.9.2-1
 module load netCDF-Fortran/4.6.1
 # module load CMake         # not needed in rundales
 # module load HDF5/1.14.2   # implicitly loaded
 # module load netCDF/4.9.2  # implicitly loaded



export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun --cpu-bind=map_cpu:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 ./dales4 > runlog.txt
###srun ./dales4 > stdlog.txt
###python3 ./change_nml.py
