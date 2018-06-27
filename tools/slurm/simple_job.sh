#!/bin/bash

#SBATCH --output=out
#SBATCH --error=err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1

# TODO: add header and description

mpirun -np 1 ./bin/fluid
