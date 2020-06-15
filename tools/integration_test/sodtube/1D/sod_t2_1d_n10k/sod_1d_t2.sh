#!/bin/bash

#SBATCH --job-name=sod_1d_t2r
#SBATCH --output=sod_1d_t2r.out
#SBATCH --error=sod_1d_t2r.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

BUILD_PATH=../../../../../build/app/ 

mpirun -np 1 ./${BUILD_PATH}/id_generators/sodtube_1d_generator sod_test2_n10000.par 
mpirun -np 1 ./${BUILD_PATH}/drivers/hydro_1d sod_test2_n10000.par 
