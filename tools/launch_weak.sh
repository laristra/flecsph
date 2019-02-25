#!/bin/bash

# NUMBER OF PARTICLES SQUARED
START_NPART=150
# OMP THREADS PER MPI PROCESS
OMP_THREADS=18
CORES=18
# MPI PROCESSES PER NODES
# DARWIN = 2 CPU * 18 CORES per NODES = 2 THREADS PER CORE
PER_NODE=2
num=1
# max 6=32 7=64 8=128 9=256 10=512
max=10
NITER=100
NOUTPUT=1000

if [ -z $1 ]
then
  echo "Error. Usage: "
  echo "./launch_weak.sh [sodtube,noh,sedov]"
  exit 1
fi

# Default = sodtube
TEST=0
INTPUT=""
if [[ "sodtube" == *$1* ]]; then
  echo "Running sodtube test"
  TEST=1
  INPUT="sodtube"
elif [[ "noh" == *$1* ]]; then
  echo "Running NOH test"
  TEST=2
  INPUT="noh"
elif [[ "sedov" == *$1* ]]; then
  echo "Running Sedov test"
  TEST=3
  INPUT="sedov"
fi

if [ $TEST -eq 0 ]; then
  echo "Error. Invalid test name"
  exit 1
fi


# Parameter files for all the tests
SODTUBE="
  #
  # Sodtube test
  #
  # initial data
  poly_gamma = 1.4        \# polytropic index
  equal_mass = yes        \# determines whether equal mass particles are used or equal separation
  sph_eta = 1.2
  lattice_type = 1        \# 0:rectangular, 1:hcp, 2:fcc
  # evolution
  sph_kernel = \"Wendland C6\"
  initial_dt = 1.0
  sph_variable_h = yes
  adaptive_timestep = yes
  timestep_cfl_factor = 0.25
  initial_iteration = 0
  final_iteration = $NITER
  out_screen_every = 20
  out_scalar_every = 20
  out_h5data_every = $NOUTPUT
  out_diagnostic_every = 1
  external_force_type=\"walls:xyz\"
  zero_potential_poison_value = 0.0
  extforce_wall_steepness = 1e12
  extforce_wall_powerindex = 5"

NOH="
#
# Noh collapse, rebounce & standing shock test
#
# initial data
  poly_gamma = 1.6666667      # polytropic index
  rho_initial = 1.0
  pressure_initial = 1.0e-6
  sphere_radius = 1.0
  sph_eta = 1.2
  lattice_type = 2         # 0:rectangular, 1:hcp, 2:fcc, 3:spherical
                           # (in 2d both hcp and fcc are triangular)
  domain_type = 1          # 0:box, 1:sphere
  # box_length = 1.0
  # box_width  = 1.3
  # box_height = 0.8

# evolution parameters:
  sph_kernel = \"quintic spline\"
  initial_dt = 2.e-3  # TODO: better use Courant factor X sph_separation
  final_iteration = $NITER
  final_time = 10.0
  out_screen_every = 20
  out_scalar_every = 20
  out_h5data_every = $NOUTPUT
  sph_variable_h = yes
  adaptive_timestep = yes
  timestep_cfl_factor = 0.25"

SEDOV="
  #
  # Sedov blast wave test
  #
  # initial data

  poly_gamma = 1.4         # polytropic index
  rho_initial = 1.0
  pressure_initial = 1.0e-7
  sphere_radius = 1.0
  sph_eta = 1.2
  sedov_blast_energy = 1.0
  sedov_blast_radius = 1.0 # in units of particle separation
  lattice_type = 2         # 0:rectangular, 1:hcp, 2:fcc, 3:spherical
                           # (in 2d both hcp and fcc are triangular)
  domain_type = 1          # 0:box, 1:sphere

  # evolution parameters:
  sph_kernel = \"quintic spline\"
  initial_dt = 2.e-3  # TODO: better use Courant factor X sph_separation
  final_iteration = $NINTER
  out_screen_every = 20
  out_scalar_every = 20
  out_h5data_every = $NOUTPUT

  sph_variable_h = yes
  adaptive_timestep = yes
  timestep_cfl_factor = 0.01"

TYPE=$SODTUBE
if [ $TEST -eq 2 ]; then
  TYPE=$NOH
elif [ $TEST -eq 3 ]; then
  TYPE=$SEDOV
fi


MAX_NODES=$(bc <<< "2^($max-2)")
MAX_PROC=$(bc <<< "$MAX_NODES*$PER_NODE ")

dir="./${INPUT}_weak_$START_NPART"
rm -rf $dir
mkdir $dir
cd $dir

cat > task.sh <<EOL
#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH -N $MAX_NODES
#SBATCH -n $MAX_PROC
#SBATCH -c $CORES
#SBATCH --output=out
#SBATCH --error=err

export CLOG_ENABLE_STDLOG=1
export OMP_NUM_THREADS=$OMP_THREADS

EOL

for i in `seq 1 $max`; do

  NPARTS=$(bc <<< "scale=0; $START_NPART")

  # Writing the data file
  config="${INPUT}_n$NPARTS.par"
  cat > $config <<EOL
    $TYPE
    initial_data_prefix = "${INPUT}_$NPARTS"
    lattice_nx = $NPARTS        \# particle lattice linear dimension
    output_h5data_prefix = "ev_${INPUT}_$NPARTS"
EOL

  cat >> task.sh <<EOL

time mpirun -n 1 ../../../bin/id_generators/${INPUT}_3d_generator ./${INPUT}_n$NPARTS.par
echo "Computation: $num"
time mpirun -n $num --bind-to socket --map-by socket ../../../bin/drivers/hydro_3d ./${INTPUT}_n$NPARTS.par
echo "Done $num"

EOL
  num=$(($((2*$num))))
  START_NPART=$(($START_NPART+$START_NPART/4))
done
sbatch task.sh
ls
