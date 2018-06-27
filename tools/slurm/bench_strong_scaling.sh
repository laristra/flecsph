#!/bin/bash

#SBATCH --output=out_strong_scaling
#SBATCH --error=err_strong_scaling
#SBATCH --time=02:00:00
#SBATCH --extra-node-info="1:16:1"
# TODO: add description

echo "Running Strong Scaling"
START=200000

for i in `seq 1 16`
do
  echo "Working on $i"
  echo "time mpirun -np $i ./bin/tree_mpilegion $START"
  time mpirun -np $i ./bin/tree_mpilegion $START
done

