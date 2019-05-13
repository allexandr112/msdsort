#!/bin/sh
module purge
module load mpi/impi/2017.4.239
# mpiicc task.cpp -std=c++14 -o ppi -O3 -fopenmp
mpicc task.cpp -std=c++14 -o ppi -O3 -fopenmp
module purge

module purge
module load mpi/openmpi/3.1.3/gcc/8
mpicxx task.cpp -std=c++14 -o ppg -O3 -fopenmp
module purge
