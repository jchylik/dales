#!/bin/bash -l
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=20gb
#SBATCH --time=23:00:00
#SBATCH --account=UniKoeln

# number of nodes in $SLURM_NNODES (default: 1)
# number of tasks in $SLURM_NTASKS (default: 1)
# number of tasks per node in $SLURM_NTASKS_PER_NODE (default: 1)
# number of threads per task in $SLURM_CPUS_PER_TASK (default: 1)


#module avail
  module load netcdf/4.9.0_gnu  gnu/7.4.0   # !! these versions agree  (2023)
  module load openmpi
  module load cmake



#echo $PATH
export RUNLOG='./runlog.txt'
echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/rrzk/lib/netcdf/4.1.3/lib"
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/rrzk/lib/netcdf/4.1.3-gcc-4.8.2/lib"

srun -n $SLURM_NTASKS ./dales4  >${RUNLOG}



