#!/bin/bash -l
#SBATCH --partition=devel-rh7
#SBATCH --ntasks=8
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem=40gb     # --mem-per-cpu=20gb
#SBATCH --time=1:00:00
#SBATCH --account=UniKoeln

# number of nodes in $SLURM_NNODES (default: 1)
# number of tasks in $SLURM_NTASKS (default: 1)
# number of tasks per node in $SLURM_NTASKS_PER_NODE (default: 1)
# number of threads per task in $SLURM_CPUS_PER_TASK (default: 1)


# load  modules
  module load netcdf/4.9.0_gnu  gnu/7.4.0   # !! these versions agree  (2023)
  module load openmpi
  module load cmake


echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/rrzk/lib/netcdf/f_4.6.0_gnu_7.4.0/lib/"

srun -n $SLURM_NTASKS dales4



