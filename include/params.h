/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file params.h
 * @author Oleg Korobkin
 * @date June 2018
 * @brief Handle global parameters: parsing parameter files, default vals etc.
 *
 * Parameters are essential for maintaining data in scientific numerical
 * experiments. They allow to launch and keep track of a battery of experiments
 * without recompiling the code. State-of-the-art simulation codes contain
 * hundreds, sometimes thousands of parameters. Clearly, it is impossible to
 * remember what are their optimal values, which also may vary case-tocase.
 * This is why it is important to store parameters in separate files.
 *
 * For FleCSPH, we picked a simple ASCII format, which is easy to read in any
 * system. The format contains pairs of "PARAM = VALUE". There are four types
 * of parameters: integer, real, string and boolean. Integer and real
 * parameters are intialized as usual, e.g.:
 *
 *   nparticles = 250
 *   poly_gamma = 1.4
 *
 * Boolean parameters are set with "yes" or "no", and string ones are set with
 * strings, which can be enclosed in single or double quotes:
 *
 *   output_conserved = yes                # produce scalar output
 *   initial_data_prefix = "sodtube_n10k"  # name of the input file
 *
 * Comments can be added to parameter files: everything which starts with '#',
 * is considered a comment and ignored. Blank spaces and empty lines are
 * ignored.
 *
 * All parameters in a parfile have exactly the same name as in the code, to
 * avoid confusion. Parameters are read-only (const references) in the param::
 * namespace. It is also possible to #define a parameter instead for optimized
 * performance -- in this case, however, this parameter needs to be commented out
 * in the parameter file.
 *
 * To introduce a new parameter:
 *  - add its declaration below using DECLARE_PARAM or DECLARE_STRING_PARAM
 *    macro (see below);
 *  - add REAL_NUMERIC_PARAM, READ_BOOLEAN_PARAM or READ_STRING_PARAM
 *    to the set_param() function -- see examples below.
 *
 * TODO: introduce a set of macros to do both declaration and reading:
 *       REGISTER_BOOLEAN_PARAM, REGISTER_INTEGER_PARAM,
 *       REGISTER_STRING_PARAM, REGISTER_REAL_PARAM
 *       instead of DECLARE/READ pair
 */

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <cstdbool>
#include "cinchlog.h"
#include "mpi.h"

#ifndef PARAMS_H
#define PARAMS_H

//////////////////////////////////////////////////////////////////////
#define DECLARE_PARAM(PTYPE,PNAME,PDEF) \
  PTYPE _##PNAME = (PDEF); \
  const PTYPE & PNAME = _##PNAME;

#define STRING_MAXLEN 511
#define DECLARE_STRING_PARAM(PNAME,PDEF) \
  char _##PNAME[STRING_MAXLEN] = (PDEF); \
  const char * PNAME = _##PNAME;

#define DO_QUOTE(str) #str
#define QUOTE(str) DO_QUOTE(str)

#define READ_BOOLEAN_PARAM(PNAME) \
  if (param_name == QUOTE(PNAME)) { \
    (_##PNAME) = lparam_value; unknown_param=false;}

#define READ_NUMERIC_PARAM(PNAME) \
  if (param_name == QUOTE(PNAME)) { \
    iss >> (_##PNAME); unknown_param=false;}

#define READ_STRING_PARAM(PNAME) \
  if (param_name == QUOTE(PNAME)) { \
    strcpy(_##PNAME, str_value.c_str()); unknown_param = false;}

#define SET_PARAM(PNAME,EXPR) _##PNAME = (EXPR);
// TODO: the macro above won't work for strings!

//////////////////////////////////////////////////////////////////////

namespace param {

//
// Parameters controlling timestepping and iterations
//
//- initial iteration
#ifndef initial_iteration
  DECLARE_PARAM(int64_t,initial_iteration,0)
#endif

//- final iteration (= total iterations + 1, if counting from 0)
#ifndef final_iteration
  DECLARE_PARAM(int64_t,final_iteration,10)
#endif

#ifndef initial_time
  DECLARE_PARAM(int64_t,initial_time,0)
#endif

#ifndef final_time
  DECLARE_PARAM(int64_t,final_time,1.0)
#endif

//- inital timestep
#ifndef initial_dt
  DECLARE_PARAM(double,initial_dt,0.001)
#endif

//
// Parameters related to particle number and density
//
//- global number of particles
#ifndef nparticles
  DECLARE_PARAM(int64_t,nparticles,1000)
#endif

//- particle lattice linear dimension
#ifndef lattice_nx
  DECLARE_PARAM(int64_t,lattice_nx,100)
#endif

//- SPH eta parameter, eta = h (rho/m)^1/D (Rosswog'09, eq.51)
#ifndef sph_eta
  DECLARE_PARAM(double,sph_eta,1.5)
#endif

//- if smoothing length is constant: value of the smoothing length
#ifndef sph_smoothing_length
  DECLARE_PARAM(double,sph_smoothing_length,-1.0) // POISONED DEFAULT
#endif

//- minimum interparticle distance (for some initial data problems)
#ifndef sph_separation
  DECLARE_PARAM(double,sph_separation,-1.0) // POISONED DEFAULT
#endif

//- which kernel to use
#ifndef sph_kernel
  DECLARE_STRING_PARAM(sph_kernel,"Wendland C2")
#endif

//- sinc kernel power index
#ifndef sph_sinc_index
  DECLARE_PARAM(double,sph_sinc_index,4.0)
#endif

//- if true, recompute (uniform) smoothing length every timestep 
//  h = average { sph_eta (m/rho)^1/D } (Rosswog'09, eq.51)
# ifndef sph_update_uniform_h
  DECLARE_PARAM(bool, sph_update_uniform_h,false)
# endif 

//
// Geometric parameters
//
//- rectangular configuration parameters (e.g. sodtube in 2D/3D)
#ifndef box_length
  DECLARE_PARAM(double,box_length,3.0)
#endif

#ifndef box_width
  DECLARE_PARAM(double,box_width,1.0)
#endif

#ifndef box_height
  DECLARE_PARAM(double,box_height,1.0)
#endif

#ifndef sphere_radius
  DECLARE_PARAM(double,sphere_radius,1.0)
#endif

//
// Boundary conditions
//
//- TODO: add description
#ifndef do_boundaries
  DECLARE_PARAM(bool,do_boundaries,false)
#endif

//- TODO: add description
#ifndef stop_boundaries
  DECLARE_PARAM(bool,stop_boundaries,false)
#endif

//- TODO: add description
#ifndef reflect_boundaries
  DECLARE_PARAM(bool,reflect_boundaries,false)
#endif

//
// I/O parameters
//
//- file prefix for input and intiial data file[s]
#ifndef initial_data_prefix
  DECLARE_STRING_PARAM(initial_data_prefix,"initial_data")
#endif

//- file prefix for HDF5 output data file[s]
#ifndef output_h5data_prefix
  DECLARE_STRING_PARAM(output_h5data_prefix,"output_data")
#endif

//- screen output frequency
#ifndef out_screen_every
  DECLARE_PARAM(int32_t,out_screen_every,1)
#endif

//- scalar reductions output frequency
#ifndef out_scalar_every
  DECLARE_PARAM(int32_t,out_scalar_every,10)
#endif

//- HDF5 output frequency
#ifndef out_h5data_every
  DECLARE_PARAM(int32_t,out_h5data_every,10)
#endif

//- produce separate HDF5 file per iteration
#ifndef out_h5data_separate_iterations
  DECLARE_PARAM(bool,out_h5data_separate_iterations,false)
#endif

//
// Viscosity and equation of state
//
//- which equation of state to use?
//  * "ideal fluid" (default)
//  * "polytropic"
//  * "white dwarf"
#ifndef eos_type
  DECLARE_STRING_PARAM(eos_type,"ideal fluid")
#endif

//- polytropic index
#ifndef poly_gamma
  DECLARE_PARAM(double,poly_gamma,1.4)
#endif

//- artificial viscosity: parameter alpha (Rosswog'09, eq.59)
#ifndef sph_viscosity_alpha
  DECLARE_PARAM(double,sph_viscosity_alpha,1.0)
#endif

//- artificial viscosity: parameter beta
#ifndef sph_viscosity_beta
  DECLARE_PARAM(double,sph_viscosity_beta,2.0)
#endif

//- artificial viscosity: parameter eta
#ifndef sph_viscosity_epsilon
  DECLARE_PARAM(double,sph_viscosity_epsilon,0.01)
#endif

//
// Gravity-related parameters
//
//- mac'n'cheese acceptance criteria
# ifndef fmm_macangle
  DECLARE_PARAM(double,fmm_macangle,0.0)
# endif

//- maximum mass per cell
# ifndef fmm_max_cell_mass
  DECLARE_PARAM(double,fmm_max_cell_mass, 1.0e-4)
# endif

//
// Drag force parameters
//
// HL : These parameters are used to relax star from 
//      initial star in both single and binary system.
//      Drag force is applied to acceleration computation
//      during beginning of steps.
//      If we have some IDs that do not require relaxation,
//      this can be neglected
//
//- Do drag froce boolean
#ifndef do_drag
  DECLARE_PARAM(bool,do_drag,false)
#endif

//- relaxation steps
# ifndef relax_steps
  DECLARE_PARAM(int,relax_steps,10)
# endif

//- Drag force coefficients. 
# ifndef drag_coeff
  DECLARE_PARAM(double,drag_coeff,1.e-6)
# endif


//
// Parameters for external acceleration
//
//- which external force to apply?
//  * "none" (default)
#ifndef external_force_type
  DECLARE_STRING_PARAM(external_force_type,"none")
#endif

# ifndef zero_potential_poison_value
  DECLARE_PARAM(double,zero_potential_poison_value, 0.0)
# endif

# ifndef extforce_wall_powerindex
  DECLARE_PARAM(double,extforce_wall_powerindex, 5.0)
# endif

# ifndef extforce_wall_steepness
  DECLARE_PARAM(double,extforce_wall_steepness, 1e6)
# endif

# ifndef thermokinetic_formulation
  DECLARE_PARAM(bool,thermokinetic_formulation, true)
# endif

//
// Specific apps
//
/// number of Sodtest to run (1..5)
#ifndef sodtest_num
  DECLARE_PARAM(unsigned short,sodtest_num,1)
#endif

// equal mass or equal particle separation switch for sodtube
#ifndef equal_mass
  DECLARE_PARAM(bool,equal_mass,true)
#endif

// characteristic density for an initial conditions
# ifndef rho_initial
  DECLARE_PARAM(double,rho_initial,1.0)
# endif

// characteristic pressure
# ifndef pressure_initial
  DECLARE_PARAM(double,pressure_initial,1.0)
# endif

// characteristic specific internal energy
# ifndef uint_initial
  DECLARE_PARAM(double,uint_initial,1.0)
# endif

// in Sedov test: total injected blast enregy
# ifndef sedov_blast_energy
  DECLARE_PARAM(double,sedov_blast_energy,1.0)
# endif

// in Sedov test: radius of energy injection
// (in units of particle separation)
# ifndef sedov_blast_radius
  DECLARE_PARAM(double,sedov_blast_radius,1.0)
# endif

// in Sedov test: set the lattice types
# ifndef lattice_type
  DECLARE_PARAM(int,lattice_type,0)
# endif

// in Sedov test: set the lattice types
# ifndef domain_type
  DECLARE_PARAM(int,domain_type,0)
# endif

// in several tests: initial velocity of the flow
# ifndef flow_velocity
  DECLARE_PARAM(double,flow_velocity,0.0)
# endif

//
// Airfoil parameters
//
# ifndef airfoil_size
  DECLARE_PARAM(double,airfoil_size, 2.0)
# endif

# ifndef airfoil_thickness
  DECLARE_PARAM(double,airfoil_thickness, 0.05)
# endif

# ifndef airfoil_camber
  DECLARE_PARAM(double,airfoil_camber, 0.1)
# endif

# ifndef airfoil_anchor_x
  DECLARE_PARAM(double,airfoil_anchor_x, -1.0)
# endif

# ifndef airfoil_anchor_y
  DECLARE_PARAM(double,airfoil_anchor_y, 0.0)
# endif

# ifndef airfoil_attack_angle
  DECLARE_PARAM(double,airfoil_attack_angle, 0.0)
# endif

// ---

/*!
   trim whitespace from the begging and end of line
 */
std::string trim(const std::string& str) {
  using namespace std;
  const size_t strBegin = str.find_first_not_of(" \t");
  if (strBegin == string::npos)
    return "";
  const size_t strEnd = str.find_last_not_of(" \t");
  return str.substr(strBegin, strEnd - strBegin + 1);
}


/**
 * @brief Sets a global parameter param_name to the value, given by its
 *        string representation param_value
 */
void set_param(const std::string& param_name,
               const std::string& param_value) {
  
  // RANK/SIZE for CLOG output 
  int rank = 0;
  int size = 1; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  using namespace std;
  bool unknown_param = true;

  // strip trailing comments
  size_t i2 = param_value.find("#");
  string str_value;
  if (i2 != string::npos)
    str_value = trim(param_value.substr(0,i2));
  else
    str_value = param_value;

  // remove (single or double) quotes
  if (str_value[0]=='\'' or str_value[0]=='\"')
    str_value = str_value.substr(1,str_value.length()-2);
  istringstream iss(str_value);

  // for boolean parameters
  bool lparam_value = (str_value == "yes"
                    or str_value == "'yes'"
                    or str_value == "\"yes\""
                    or str_value == "true"
                    or str_value == "'true'"
                    or str_value == "\"true\"");
  // timestepping and iterations --------------------------------------------
# ifndef initial_iteration
  READ_NUMERIC_PARAM(initial_iteration)
# endif

# ifndef final_iteration
  READ_NUMERIC_PARAM(final_iteration)
# endif

# ifndef initial_time
  READ_NUMERIC_PARAM(initial_time)
# endif

# ifndef final_time
  READ_NUMERIC_PARAM(final_time)
# endif

# ifndef initial_dt
  READ_NUMERIC_PARAM(initial_dt)
# endif


  // particle number and density --------------------------------------------
# ifndef nparticles
  READ_NUMERIC_PARAM(nparticles)
# endif

# ifndef lattice_nx
  READ_NUMERIC_PARAM(lattice_nx)
# endif

# ifndef sph_eta
  READ_NUMERIC_PARAM(sph_eta)
# endif

# ifndef sph_smoothing_length
  READ_NUMERIC_PARAM(sph_smoothing_length)
# endif

# ifndef sph_separation
  READ_NUMERIC_PARAM(sph_separation)
# endif

# ifndef initial_data_prefix
  READ_STRING_PARAM(sph_kernel)
# endif

# ifndef sph_sinc_index
  READ_NUMERIC_PARAM(sph_sinc_index)
# endif

# ifndef sph_update_uniform_h
  READ_BOOLEAN_PARAM(sph_update_uniform_h)
# endif 

  // geometric configuration  -----------------------------------------------
# ifndef box_length
  READ_NUMERIC_PARAM(box_length)
# endif

# ifndef box_width
  READ_NUMERIC_PARAM(box_width)
# endif

# ifndef box_height
  READ_NUMERIC_PARAM(box_height)
# endif

# ifndef sphere_radius
  READ_NUMERIC_PARAM(sphere_radius)
# endif

  // boundary conditions  ---------------------------------------------------
# ifndef do_boundaries
  READ_BOOLEAN_PARAM(do_boundaries)
# endif

# ifndef stop_boundaries
  READ_BOOLEAN_PARAM(stop_boundaries)
# endif

# ifndef reflect_boundaries
  READ_BOOLEAN_PARAM(reflect_boundaries)
# endif

  // i/o parameters  --------------------------------------------------------
# ifndef initial_data_prefix
  READ_STRING_PARAM(initial_data_prefix)
# endif

# ifndef output_h5data_prefix
  READ_STRING_PARAM(output_h5data_prefix)
# endif

# ifndef out_screen_every
  READ_NUMERIC_PARAM(out_screen_every)
# endif

# ifndef out_scalar_every
  READ_NUMERIC_PARAM(out_scalar_every)
# endif

# ifndef out_h5data_every
  READ_NUMERIC_PARAM(out_h5data_every)
# endif

# ifndef out_h5data_separate_iterations
  READ_BOOLEAN_PARAM(out_h5data_separate_iterations)
# endif

  // viscosity and equation of state ----------------------------------------
# ifndef eos_type
  READ_STRING_PARAM(eos_type)
# endif

# ifndef poly_gamma
  READ_NUMERIC_PARAM(poly_gamma)
# endif

# ifndef sph_viscosity_alpha
  READ_NUMERIC_PARAM(sph_viscosity_alpha)
# endif

# ifndef sph_viscosity_beta
  READ_NUMERIC_PARAM(sph_viscosity_beta)
# endif

# ifndef sph_viscosity_epsilon
  READ_NUMERIC_PARAM(sph_viscosity_epsilon)
# endif

  // gravity-related  -------------------------------------------------------
# ifndef fmm_macangle
  READ_NUMERIC_PARAM(fmm_macangle)
# endif

# ifndef fmm_max_cell_mass
  READ_NUMERIC_PARAM(fmm_max_cell_mass)
# endif

// Drag force parameters ---------------------------------------------------
#ifndef do_drag
  READ_BOOLEAN_PARAM(do_drag)
#endif

# ifndef relax_steps
  READ_NUMERIC_PARAM(relax_steps)
# endif

# ifndef drag_coeff
  READ_NUMERIC_PARAM(drag_coeff)
# endif

// external force  --------------------------------------------------------
#ifndef external_force_type
  READ_STRING_PARAM(external_force_type)
#endif

# ifndef zero_potential_poison_value
  READ_NUMERIC_PARAM(zero_potential_poison_value)
# endif

# ifndef extforce_wall_powerindex
  READ_NUMERIC_PARAM(extforce_wall_powerindex)
# endif

# ifndef extforce_wall_steepness
  READ_NUMERIC_PARAM(extforce_wall_steepness)
# endif

# ifndef thermokinetic_formulation
  READ_BOOLEAN_PARAM(thermokinetic_formulation)
# endif

  // specific apps  ---------------------------------------------------------
# ifndef sodtest_num
  READ_NUMERIC_PARAM(sodtest_num)
# endif

#ifndef equal_mass
  READ_BOOLEAN_PARAM(equal_mass)
#endif

# ifndef rho_initial
  READ_NUMERIC_PARAM(rho_initial)
# endif

# ifndef pressure_initial
  READ_NUMERIC_PARAM(pressure_initial)
# endif

# ifndef uint_initial
  READ_NUMERIC_PARAM(uint_initial)
# endif

# ifndef sedov_blast_energy
  READ_NUMERIC_PARAM(sedov_blast_energy)
# endif

# ifndef sedov_blast_radius
  READ_NUMERIC_PARAM(sedov_blast_radius)
# endif

# ifndef lattice_type
  READ_NUMERIC_PARAM(lattice_type)
# endif

# ifndef domain_type
  READ_NUMERIC_PARAM(domain_type)
# endif

# ifndef flow_velocity
  READ_NUMERIC_PARAM(flow_velocity)
# endif

  // airfoil parameters  ----------------------------------------------------
# ifndef airfoil_size
  READ_NUMERIC_PARAM(airfoil_size)
# endif

# ifndef airfoil_thickness
  READ_NUMERIC_PARAM(airfoil_thickness)
# endif

# ifndef airfoil_camber
  READ_NUMERIC_PARAM(airfoil_camber)
# endif

# ifndef airfoil_anchor_x
  READ_NUMERIC_PARAM(airfoil_anchor_x)
# endif

# ifndef airfoil_anchor_y
  READ_NUMERIC_PARAM(airfoil_anchor_y)
# endif

# ifndef airfoil_attack_angle
  READ_NUMERIC_PARAM(airfoil_attack_angle)
# endif

  // unknown parameter -------------------------------
  if (unknown_param) {
    clog(error) << "ERROR: unknown parameter " << param_name << endl;
    exit(2);
  }

  rank || clog(trace) << param_name << ": " << param_value << endl;
}

/**
 * @brief Parameter file parser
 *
 * Simulation parameter file contains strings of pairs:
 *  param = value
 *
 * assigining values to specific parameters. Parameters can have four
 * different types: logical (boolean), integer, real and string.
 * Parameter file can contain comments. Comments begin with '#';
 * everything after '#' is ignored. Strings can be embraced in single
 * quotes (') or double quotes ("). Logical parameters are initialized
 * with "yes" or "no".
 *
 * Example parameter file:
 * --->>> sodtube_np10k.par ---------------
  # Sodtube test #1 for 10000 particles
  nparticles= 10000   # global number of particles
  poly_gamma = 1.4    # polytropic index
  sodtest_num  = 1    # which Sod test to generate
  output_singlefile = yes
  initial_data_prefix  = "sodtube_np10k"
 * ---<<<  --------------------------------
 */
void read_params(const char * parfile) {
  using namespace std;
  ifstream infile;
  string line;

  // attempt to open the file
  infile.open (parfile);
  if (!infile) {
    cerr << "ERROR: Unable to open parameter file 'input.par'" << endl;
    exit(1);
  }


  for (int ln=1; std::getline(infile,line); ++ln) {

    // skip comments (lines starting with '#' at any position)
    // and blank lines
    bool is_blank = true, is_comment = false;
    for (size_t i=0;i<line.length();i++) {
      char c = line[i];
      is_comment = (c == '#');
      is_blank = (c == ' ' || c == '\t');
      if (is_comment or not is_blank) break;
    }
    if (is_comment or is_blank) continue;

    // line must be in the form:
    // lhs = rhs
    size_t eqpos = line.find("=");
    if (eqpos == string::npos) {
      cerr << "ERROR in parameter file " << parfile
           << ", line #" << ln << endl;
      exit(1);
    }
    string lhs = line.substr(0,eqpos);
    string rhs = line.substr(eqpos+1);

    // strip whitespace
    lhs = trim(lhs);
    rhs = trim(rhs);
    if (lhs.length() == 0 or rhs.length() == 0) {
      cerr << "ERROR in parameter file " << parfile
           << ", line #" << ln << ": wrong format" << endl;
      exit(1);
    }
    set_param(lhs,rhs);
  }
  infile.close();

}

/**
 * @brief MPI parameter file reader
 * @todo  Use FleCSI infrastructure instead (e.g. Flecsi_Sim_IO?)
 */
void mpi_read_params(const char * parameter_file) {
  const int MAXLEN = 2048;
  char buffer[MAXLEN];
  char * parfile = buffer;
  int len, rank, size, parfile_free = 0;

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Status status;

  // the char* array pointer is only valid on rank 0;
  // broadcast parfile name from rank 0 to other ranks
  if (rank == 0) {
    strcpy(parfile,parameter_file);
    len = strlen(parfile);
  }
  MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(parfile, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);

  rank || clog(trace) << "Parameter file name on rank " << rank << " over "<<
              size << ": " << parfile << std::endl << std::flush;

  // queue ranks to read the parfile sequentially;
  // wait for a message from previous rank, unless this is rank 0
  if (rank > 0)
    MPI_Recv(&parfile_free,1,MPI_INT,rank-1,1,MPI_COMM_WORLD,&status);

  // read parameters
  read_params (parfile);

  // inform the next rank that the file is free to read
  if (rank + 1 < size)
    MPI_Send(&parfile_free,1,MPI_INT,rank+1,1,MPI_COMM_WORLD);

} // mpi_read_param

} // namespace params

#endif // PARAMS_H
