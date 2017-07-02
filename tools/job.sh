#!/bin/bash

#SBATCH -n 32
#SBATCH -t 00:30:00
#SBATCH -J flecsph

#SBATCH --output=./out
#SBATCH --error=./err

srun -n 32 perf record -g -o toto_$$.data ./bin/tree_mpilegion ../data/data_binary_8338.txt
