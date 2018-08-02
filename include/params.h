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

//- square root of the total number of particles (for 2D setups)
#ifndef sqrt_nparticles
  DECLARE_PARAM(int64_t,sqrt_nparticles,100)
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
#ifndef initial_data_prefix
  DECLARE_STRING_PARAM(sph_kernel,"Wendland quintic")
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

//
// Viscosity and equation of state
//
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
// Specific apps
//
/// number of Sodtest to run (1..5)
#ifndef sodtest_num
  DECLARE_PARAM(unsigned short,sodtest_num,1)
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
                    or str_value == "\"yes\"");
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

# ifndef sqrt_nparticles
  READ_NUMERIC_PARAM(sqrt_nparticles)
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

  // viscosity and equation of state ----------------------------------------
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

  // specific apps  ---------------------------------------------------------
# ifndef sodtest_num
  READ_NUMERIC_PARAM(sodtest_num)
# endif

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
  // unknown parameter -------------------------------
  if (unknown_param) {
    clog_one(fatal) << "ERROR: unknown parameter " << param_name << endl;
    exit(2);
  }

  clog_one(trace) << param_name << ": " << param_value << endl;
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

  clog(trace) << "Parameter file name on rank " << rank
              << ": " << parfile << std::endl << std::flush;

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
