#!/bin/bash

#SBATCH -t 10:00:00

#SBATCH -N 4
#SBATCH -n 8
#SBATCH -c 8

#SBATCH --output=out
#SBATCH --error=err

#SBATCH --gres=gpu:2 
# TODO: add header and description

source ~/flecsi/load_flecsi_modules.sh
#module load maqao 

time mpirun -np 1 ./sedov_generator 

for i in `seq 1 8`;
do

  export OMP_NUM_THREADS=1
  /usr/bin/time -f "TIME FOR $i/1 %e" mpirun -np $i ./sedov 

  #export OMP_NUM_THREADS=2
  #/usr/bin/time -f "TIME FOR $i/2 %e" mpirun -np $i ./sedov 

  #export OMP_NUM_THREADS=4
  #/usr/bin/time -f "TIME FOR $i/4 %e" mpirun -np $i ./sedov 

  #export OMP_NUM_THREADS=8
  #/usr/bin/time -f "TIME FOR $i/8 %e" mpirun -np $i ./sedov 

done 

time mpirun -np 1 ./sedov_generator 200 
for i in `seq 1 8`;
do
  export OMP_NUM_THREADS=1
  /usr/bin/time -f "TIME FOR $i/1 %e" mpirun -np $i ./sedov 
done 

time mpirun -np 1 ./sedov_generator 300 
for i in `seq 1 8`;
do
  export OMP_NUM_THREADS=1
  /usr/bin/time -f "TIME FOR $i/1 %e" mpirun -np $i ./sedov 
done 

