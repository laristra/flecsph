#!/bin/bash
# NUMBER OF PARTICLES :: calculated based on the dimension of the test. (ONLY DIM=3)
#                     :: currently aiming for 20,000 particles per core
#                     :: int(cbrt(20000*cores)+1)
NPART_PER_CORE=20000
CORE_MAX=36
# max_ranks 0=1 1=2 2=4 3=8 :: max number of ranks
max_ranks=5
# max_omp 0=1 1=2 2=4 3=8 :: max number of omp_threads
max_omp=5
# number of iterations to perform evolution over
NITER=100
# iteration number for writing to file :: this needs to larger than NITER in order to not time read/write
NOUTPUT=1000

if [ -z $1 ]
then
  echo "Error. Usage: "
  echo "./launch_strong_openMP.sh [sodtube,noh,sedov]"
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
  poly_gamma   = 1.4      # polytropic index
  sodtest_num  = 1        # which test to generate (1..6)
  equal_mass   = yes      # determines whether equal mass particles are used or equal separation
  lattice_type = 0        # 0:rectangular, 1:hcp, 2:fcc
  domain_type  = 0        # 0:box
  box_length   = 3.0
  box_width    = 1.0
  box_height   = 1.0
  # evolution
  sph_kernel          = \"Wendland C6\"
  sph_eta             = 1.5
  initial_dt          = 1.e-10
  sph_variable_h      = yes
  adaptive_timestep   = yes
  timestep_cfl_factor = 0.1

  final_iteration       = $NITER
  out_screen_every      = 20
  out_scalar_every      = 20
  out_h5data_every      = $NOUTPUT
  out_diagnostic_every  = 1
  external_force_type   = \"walls:xyz\"
  zero_potential_poison_value = 0.0
  extforce_wall_steepness     = 1e10
  extforce_wall_powerindex    = 3"

NOH="
#
# Noh collapse, rebounce & standing shock test
#

# geometry
  domain_type  = 0         # 0:box, 1:sphere
  lattice_type = 0         # 0:rectangular, 1:hcp, 2:fcc
  box_length   = 1.0
  box_width    = 1.0
  box_height   = 1.0
# EoS
  eos_type   = \"ideal fluid\"            # EoS: ideal fluid, polytropic, white dwarf, ...
  poly_gamma = 1.66667                    # Polytropic exponent
# density and pressure
  rho_initial      = 1.0                  # Initial density
  pressure_initial = 1.e-7                # Initial pressure
# Noh
  noh_infall_velocity = 0.1               # Infall velocity in Noh test
# evolution parameters:
  sph_kernel          = \"Wendland C6\"   # Name of SPH kernel, e.g. Wendland
  sph_eta             = 1.6               # Smoothing length, e.g. 1.6
  sph_variable_h      = yes               # Adaptive smoothing length switch
  initial_dt          = 1.e-10            # Initial timestep
  timestep_cfl_factor = 0.1               # Courant–Friedrichs–Lewy factor
  adaptive_timestep   = yes               # Adaptive timestep switch
  thermokinetic_formulation = yes         # Evolves total (=thermal+potential+kinetic) energy
# output
  final_iteration = $NITER
  final_time = 10.0
  out_screen_every = 20
  out_scalar_every = 20
  out_h5data_every = $NOUTPUT"

SEDOV="
 #
 # Sedov blast wave test
 #
 # initial data
  poly_gamma       = 1.4    # polytropic index
  rho_initial      = 1.0
  pressure_initial = 1.0e-7
  sph_eta          = 1.2
  sedov_blast_energy = 1.e-6
  sedov_blast_radius = 0.06               # in units of particle separation
  lattice_type = 0                        # 0:rectangular, 1:hcp, 2:fcc, 3:spherical
  domain_type  = 0                        # 0:box, 1:sphere
  box_length   = 2.0
  box_width    = 2.0
  box_height   = 2.0
 # evolution parameters:
  sph_kernel = \"Wendland C6\"
  initial_dt          = 1.e-10            # Initial timestep
  timestep_cfl_factor = 0.1               # Courant–Friedrichs–Lewy factor
  adaptive_timestep   = yes               # Adaptive timestep switch
  final_iteration     = $NITER
  out_screen_every    = 20
  out_scalar_every    = 20
  out_h5data_every    = $NOUTPUT
  sph_variable_h      = yes"

TYPE=$SODTUBE
if [ $TEST -eq 2 ]; then
  TYPE=$NOH
elif [ $TEST -eq 3 ]; then
  TYPE=$SEDOV
fi

NPART_MAX=$(bc -l <<< "e(l($NPART_PER_CORE*$CORE_MAX)/3)+1.5")
NPART_MAX=$(bc <<< "$NPART_MAX/1")
dir="./${INPUT}_strong_openMP_$NPART_MAX"
rm -rf $dir
mkdir $dir
cd $dir

# Writing the data file
config="${INPUT}_n${NPART_MAX}.par"
cat > $config <<EOL
  $TYPE
  initial_data_prefix  = "${INPUT}_${NPART_MAX}"
  lattice_nx           = ${NPART_MAX}               # particle lattice linear dimension
  output_h5data_prefix = "ev_${INPUT}_${NPART_MAX}"
EOL

cat > task.sh <<EOL
#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 18
#SBATCH --output=out
#SBATCH --error=err
export CLOG_ENABLE_STDLOG=1

export OMP_NUM_THREADS=$CORE_MAX
time mpirun -n 1 ${PATH_TO_FLECSPH_DIR}/build/app/id_generators/${INPUT}_3d_generator ./${INPUT}_n${NPART_MAX}.par
EOL

for i in `seq 0 $max_ranks`; do
  for j in `seq 0 $max_omp`; do
    num_ranks=$(bc <<< "scale=0; 2^$i")
    num_omp=$(bc <<< "scale=0; 2^$j")
    num_cores=$(bc <<< "scale=0; $num_ranks*$num_omp")
    if [ $num_cores -le 36 ]; then
  cat >> task.sh <<EOL

export OMP_NUM_THREADS=$num_omp
echo "Computation: RANKS ${num_ranks}, OMP_THREADS $num_omp"
time mpirun -n $num_ranks --bind-to socket --map-by socket ${PATH_TO_FLECSPH_DIR}/build/app/drivers/hydro_3d ./${INPUT}_n${NPART_MAX}.par
echo "Done: RANKS ${num_ranks}, OMP_THREADS $num_omp"
rm ev_${INPUT}_${NPART_MAX}.h5part

EOL
    fi
  done
done
#sbatch task.sh
#ls
