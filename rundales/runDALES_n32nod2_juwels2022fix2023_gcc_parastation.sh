#!/bin/bash -x
#SBATCH --account=rcongm
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --output=./slurm_dales-%j.out
#SBATCH --error=./slurm_error-%j.out
#SBATCH --time=4:00:00
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
# module use /gpfs/software/juwels/otherstages
# module load Stages/2023  #< default, not needed


module load Stages/2022

module load GCC/11.2.0  ParaStationMPI/5.5.0-1
module load netCDF-Fortran/4.5.3
module load HDF5/1.12.1    # <- loadd automatically
module load netCDF/4.8.1   # <- laoded automatically
module load CMake



# ---- start job  --------------------
echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS

# adding a busy file 
echo "runDALES job ${SLURM_JOB_ID} is busy ... don't interfere!" > $busyfile


srun ${DALES_EXECUTABLE}>$RUNLOG
#or:  srun -n $SLURM_NTASKS ${DALES_EXECUTABLE}>$RUNLOG

rm $busyfile

