#!/bin/sh
module purge
module load mpi/impi/2017.4.239
mpiicc task.cpp -o ppi -O3 -fopenmp
module purge

module purge
module load mpi/openmpi/3.1.3/gcc/8
mpicxx task.cpp -o ppg -O3 -fopenmp
module purge
