#!/bin/bash -x
#SBATCH --account=rcongm
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --output=./slurm_dales-%j.out
#SBATCH --error=./slurm_error-%j.out
#SBATCH --time=24:00:00
#SBATCH --partition=batch


# set a few variables
    export DALES_EXECUTABLE='dales4'
    export busyfile='runDALES_busy.txt'
    export RUNLOG='./runlog.txt'
    export MSGFILE='dales_mosaic_ayil.txt'

# Note: The current working directory at this point is
# the directory where sbatch was executed.

echo "Starting run"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# --- load modules ------------------

 module load Stages/2024
 module load GCC/12.3.0   ParaStationMPI/5.9.2-1
 module load netCDF-Fortran/4.6.1
 # module load CMake         # not needed in rundales
 # module load HDF5/1.14.2   # implicitly loaded
 # module load netCDF/4.9.2  # implicitly loaded



# ---- start job  --------------------
echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS
echo "runDALES job ${SLURM_JOB_ID} is busy ... don't interfere!" > $busyfile


srun ${DALES_EXECUTABLE}>$RUNLOG
#or:  srun -n $SLURM_NTASKS ${DALES_EXECUTABLE}>$RUNLOG

rm $busyfile     # cleanup 

