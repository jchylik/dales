#!/bin/bash -x
#SBATCH --account=bb1339
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --output=./slurm_dales-%j.out
#SBATCH --error=./slurm_error-%j.out
#SBATCH --time=8:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL

# differences from Juwels:
# --account=rcongm      -->   --account=bb1339
# --partition=batch     -->   --partition=compute  
# --time=8:00:00        : time limit pn Levante is 8 hours
# --mail-type=FAIL      : notify us if crash



  # set a few variables
   export DALES_EXECUTABLE='dales4'
   export busyfile='runDALES_busy.txt'
   export RUNLOG='./runlog.txt'
   export MSGFILE='dales_mosaic_ayil.txt'

  # a bit of outputs 
   echo "Starting run"
   export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# --- load modules ------------------

  # all different from Juwels
   module load  openmpi/4.1.2-gcc-11.2.0
   module load  netcdf-fortran/4.5.3-openmpi-4.1.2-gcc-11.2.0 
   module load  netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0
   module load  hdf5/1.12.1-openmpi-4.1.2-gcc-11.2.0  
   module load  gcc/11.2.0-gcc-11.2.0

# ---- start job  --------------------
   echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS

  # adding a busy file 
   echo "runDALES job ${SLURM_JOB_ID} is busy ... don't interfere!" > $busyfile

  # srun command  
   srun ${DALES_EXECUTABLE}>$RUNLOG
   #or:  srun -n $SLURM_NTASKS ${DALES_EXECUTABLE}>$RUNLOG

  # end, remove busy file 
  rm $busyfile

