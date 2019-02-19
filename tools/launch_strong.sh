#!/bin/bash 

# NUMBER OF PARTICLES SQUARED 
START_NPART=500

# OMP THREADS PER MPI PROCESS
OMP_THREADS=18
CORES=18

# MPI PROCESSES PER NODES 
# DARWIN = 2 CPU * 18 CORES per NODES = 2 THREADS PER CORE 
PER_NODE=2

num=1
# max 6=32 7=64 8=128 9=256 10=512
max=7

MAX_NODES=32
MAX_PROC=64

NITER=50
NOUTPUT=1000

for i in `seq 1 $max`; do
  pwd
  dir_mpi="mpi${num}_${OMP_THREADS}_strong"
  echo "Directory: $dir_mpi"
  rm -rf $dir_mpi
  mkdir $dir_mpi
  cd $dir_mpi
  pwd
  NODES=$(($num/$PER_NODE))
  if [ $NODES -eq 0 ]; then 
    NODES=1
  fi
  NPARTS=$(bc <<< "scale=0; $START_NPART")

  # Create the task 
  echo "#!/bin/bash" >> task.sh
  echo "#SBATCH --time=10:00:00" >> task.sh
  echo "#SBATCH -N $MAX_NODES" >> task.sh
  echo "#SBATCH -p scaling" >> task.sh 
  echo "#SBATCH -n $MAX_PROC" >> task.sh
  echo "#SBATCH -c $CORES" >> task.sh
#  echo "#SBATCH --nodelist=cn[314-329]" >> task.sh
  echo "#SBATCH --output=out$NPARTS" >> task.sh
  echo "#SBATCH --error=err$NPARTS" >> task.sh
  echo "export CLOG_ENABLE_STDLOG=1" >> task.sh
  echo "export OMP_NUM_THREADS=$OMP_THREADS" >> task.sh
  echo "time mpirun -n 1 ../../../bin/id_generators/sodtube_3d_generator ./sodtube_n$NPARTS.par" >> task.sh
  echo "time mpirun -n $num ../../../bin/drivers/hydro_3d ./sodtube_n$NPARTS.par" >> task.sh

  config="sodtube_n$NPARTS.par"
  # Create the config file 
  echo "#" >> $config 
  echo "# Sodtube test" >> $config 
  echo "#">> $config 
  echo "# initial data">> $config 
  echo "initial_data_prefix = \"ev_st_n$NPARTS\"">> $config 
  echo "lattice_nx = $NPARTS        \# particle lattice linear dimension">> $config 
  echo "poly_gamma = 1.4        \# polytropic index">> $config 
  echo "equal_mass = yes        \# determines whether equal mass particles are used or equal separation">> $config 
  echo "sph_eta = 1.2">> $config 
  echo "lattice_type = 1        \# 0:rectangular, 1:hcp, 2:fcc">> $config 
  echo "# evolution">> $config 
  echo "sph_kernel = \"Wendland C6\"">> $config 
  echo "initial_dt = 1.0">> $config 
  echo "sph_variable_h = yes">> $config 
  echo "adaptive_timestep = yes">> $config 
  echo "timestep_cfl_factor = 0.25">> $config 
  echo "initial_iteration = 0">> $config 
  echo "final_iteration = $NITER">> $config 
  echo "out_screen_every = 1">> $config 
  echo "out_scalar_every = 1">> $config 
  echo "out_h5data_every = $NOUTPUT">> $config 
  echo "out_diagnostic_every = 1">> $config 
  echo "output_h5data_prefix = \"ev_st_n$NPARTS\"">> $config 
  echo "external_force_type=\"walls:xyz\"">> $config
  echo "zero_potential_poison_value = 0.0">>$config
  echo "extforce_wall_steepness = 1e12">>$config
  echo "extforce_wall_powerindex = 5">>$config

  #echo -e $CONFIG >> "sodtube_n$NPARTS.par"
  #echo -e $TASK >> task.sh
  ls 
  sbatch task.sh 
  cd .. 
  num=$(($((2*$num))))
  #START_NPART=$(($START_NPART*2))
done 
