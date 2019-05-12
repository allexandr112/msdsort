# MPI MSD Sort

## How to run on Ubuntu
1. sudo apt update && sudo apt install mpich && sudo apt install libomp-dev
2. make

## Possible problems
Sometimes it can't be compiled, because it can't find path to MPI header.

In that case try to change these `#include <mpich/mpi.h>` or `#include <mpi.h>` in `task.cpp` file.