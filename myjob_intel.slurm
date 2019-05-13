#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --job-name=hpcg
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --partition=tornado
#SBATCH --verbose

. /etc/profile.d/modules.sh
module purge
module load mpi/impi/2017.4.239

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=CLOSE
export OMP_WAIT_POLICY=ACTIVE
export HWLOC_HIDE_ERRORS=1

mpiexec ./ppi 100
