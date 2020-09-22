#!/bin/bash -l
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=20gb
#SBATCH --time=120:00:00
#SBATCH --account=UniKoeln

# number of nodes in $SLURM_NNODES (default: 1)
# number of tasks in $SLURM_NTASKS (default: 1)
# number of tasks per node in $SLURM_NTASKS_PER_NODE (default: 1)
# number of threads per task in $SLURM_CPUS_PER_TASK (default: 1)


#module avail
#module load intelmpi/4.1.0 netcdf/4.1.3
 module load gnu/4.8.2 netcdf/4.1.3-gcc-4.8.2 intelmpi/2018
 module load hdf5/1.8.11 szlib
# module load cmake

#echo $PATH

echo "# of nodes, tasks:" $SLURM_NNODES $SLURM_NTASKS

#LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/rneggers/bin/netcdf-4.3.0_ifort/lib"
#LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/rrzk/lib/netcdf/4.1.3/lib"
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/rrzk/lib/netcdf/4.1.3-gcc-4.8.2/lib"

srun -n $SLURM_NTASKS ./dales4
