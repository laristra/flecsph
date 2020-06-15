#!/bin/bash 

if [[ $1 == '1D' ]]
then
  echo Running Sod shock tube 1D tests
  cd 1D/sod_t1_1d_n10k/
  sbatch sod_1d_t1.sh
  cd ../sod_t2_1d_n10k/
  sbatch sod_1d_t2.sh
  cd ../sod_t3_1d_n10k/
  sbatch sod_1d_t3.sh
  cd ../sod_t4_1d_n10k/
  sbatch sod_1d_t4.sh
  cd ../sod_t5_1d_n10k/
  sbatch sod_1d_t5.sh
elif [[ $1 == '2D' ]]
then
  echo Running Sod shock tube 2D tests
elif [[ $1 == '3D' ]]
then
  echo Running Sod shock tube 3D tests
else
  echo Running all Sod shock tube tests
fi
