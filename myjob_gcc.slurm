#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=hpcg
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --partition=tornado
#SBATCH --verbose

. /etc/profile.d/modules.sh
module purge
module load mpi/openmpi/3.1.3/gcc/8

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=CLOSE
export OMP_WAIT_POLICY=ACTIVE
export HWLOC_HIDE_ERRORS=1

mpiexec ./ppg 100000
