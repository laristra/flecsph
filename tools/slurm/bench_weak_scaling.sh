#!/bin/bash

#SBATCH --output=out_weak_scaling
#SBATCH --error=err_weak_scaling
#SBATCH --time=02:00:00
#SBATCH --extra-node-info="1:16:1"
# TODO: add header and description

echo "Running Weak Scaling"
START=200000
ADD=200000

for i in `seq 1 16`
do
  echo "Working on $i"
  echo "time mpirun -np $i ./bin/tree_mpilegion $START"
  time mpirun -np $i ./bin/tree_mpilegion $START
  ((START = START + ADD))
done

