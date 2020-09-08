/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * performance -- in this case, however, this parameter needs to be commented
 * out in the parameter file.
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

#include "log.h"
#include "mpi.h"
#include <assert.h>
#include <cstdbool>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <string>

#ifndef PARAMS_H
#define PARAMS_H
#include <boost/algorithm/string.hpp>

//////////////////////////////////////////////////////////////////////
#define DECLARE_PARAM(PTYPE, PNAME, PDEF)                                      \
  PTYPE _##PNAME = (PDEF);                                                     \
  const PTYPE & PNAME = _##PNAME;

#define STRING_MAXLEN 511
#define DECLARE_STRING_PARAM(PNAME, PDEF)                                      \
  char _##PNAME[STRING_MAXLEN] = (PDEF);                                       \
  const char * PNAME = _##PNAME;

#define DECLARE_KEYWORD_PARAM(PNAME, PDEF)                                     \
  PNAME##_keyword _##PNAME = (PDEF);                                           \
  const PNAME##_keyword & PNAME = _##PNAME;

#define DO_QUOTE(str) #str
#define QUOTE(str) DO_QUOTE(str)

#define READ_BOOLEAN_PARAM(PNAME)                                              \
  if(param_name == QUOTE(PNAME)) {                                             \
    (_##PNAME) = lparam_value;                                                 \
    unknown_param = false;                                                     \
  }

#define READ_NUMERIC_PARAM(PNAME)                                              \
  if(param_name == QUOTE(PNAME)) {                                             \
    iss >> (_##PNAME);                                                         \
    unknown_param = false;                                                     \
  }

#define READ_STRING_PARAM(PNAME)                                               \
  if(param_name == QUOTE(PNAME)) {                                             \
    strcpy(_##PNAME, str_value.c_str());                                       \
    unknown_param = false;                                                     \
  }

#define SET_PARAM(PNAME, EXPR) _##PNAME = (EXPR);
// TODO: the macro above won't work for strings!

//////////////////////////////////////////////////////////////////////
namespace param {
//
// Enums for keyword-type parameters
//

// sph_kernel keywords
typedef enum sph_kernel_keyword_enum {
  cubic_spline,
  quintic_spline,
  wendland_c2,
  wendland_c4,
  wendland_c6,
  gaussian,
  super_gaussian,
  sinc_ker
} sph_kernel_keyword;

typedef enum eos_type_keyword_enum{
  eos_ideal,
  eos_polytropic,
  eos_wd,
  eos_ppt,
  eos_no_eos
} eos_type_keyword;


// sph_viscosity keywords
typedef enum sph_viscosity_keyword_enum {
  visc_constant,
  visc_cullen
} sph_viscosity_keyword;

//////////////////////////////////////////////////////////////////////
//
// Parameters controlling timestepping and iterations
//
//- initial iteration
#ifndef initial_iteration
DECLARE_PARAM(int64_t, initial_iteration, 0)
#endif

//- final iteration (= total iterations + 1, if counting from 0)
#ifndef final_iteration
DECLARE_PARAM(int64_t, final_iteration, 0)
#endif

#ifndef initial_time
DECLARE_PARAM(int64_t, initial_time, 0)
#endif

#ifndef final_time
DECLARE_PARAM(double, final_time, 0.0)
#endif

//- inital timestep
#ifndef initial_dt
DECLARE_PARAM(double, initial_dt, 0.001)
#endif

//- timestep Courant-Friedrichs-Lewy factor (Dt/Dx)
#ifndef timestep_cfl_factor
DECLARE_PARAM(double, timestep_cfl_factor, 0.25)
#endif

//- adaptive timestepping flag
#ifndef adaptive_timestep
DECLARE_PARAM(bool, adaptive_timestep, false)
#endif

//- number of passes when computing du/dt or de/dt 
//  to accurately update the pressure (1 or 2)
#ifndef pressure_updates_number
  DECLARE_PARAM(int64_t,pressure_updates_number,1)
#endif

//
// Parameters related to particle number and density
//
//- global number of particles
#ifndef nparticles
DECLARE_PARAM(int64_t, nparticles, 1000)
#endif

//- particle lattice linear dimension
#ifndef lattice_nx
DECLARE_PARAM(int64_t, lattice_nx, 100)
#endif

//- SPH eta parameter, eta = h (rho/m)^1/D (Rosswog'09, eq.51)
#ifndef sph_eta
DECLARE_PARAM(double, sph_eta, 1.5)
#endif

//- if smoothing length is constant: value of the smoothing length
#ifndef sph_smoothing_length
DECLARE_PARAM(double, sph_smoothing_length, -1.0) // POISONED DEFAULT
#endif

//- minimum interparticle distance (for some initial data problems)
#ifndef sph_separation
DECLARE_PARAM(double, sph_separation, -1.0) // POISONED DEFAULT
#endif

//- which kernel to use
#ifndef sph_kernel
DECLARE_KEYWORD_PARAM(sph_kernel, wendland_c4)
#endif

//- sinc kernel power index
#ifndef sph_sinc_index
DECLARE_PARAM(double, sph_sinc_index, 4.0)
#endif

//- if true, recompute (uniform) smoothing length every timestep
//  h = average { sph_eta (m/rho)^1/D } (Rosswog'09, eq.51)
#ifndef sph_update_uniform_h
DECLARE_PARAM(bool, sph_update_uniform_h, false)
#endif

//- if true, the smoothing length is variable, not the same among the
// particles.
#ifndef sph_variable_h
DECLARE_PARAM(bool, sph_variable_h, false)
#endif

//
// Geometric parameters
//
//- rectangular configuration parameters (e.g. sodtube in 2D/3D)

// in various tests: sets the type of the domain (0:box, 1:sphere/circle)
#ifndef domain_type
DECLARE_PARAM(int, domain_type, 0)
#endif

#ifndef box_length
DECLARE_PARAM(double, box_length, 3.0)
#endif

#ifndef box_width
DECLARE_PARAM(double, box_width, 1.0)
#endif

#ifndef box_height
DECLARE_PARAM(double, box_height, 1.0)
#endif

#ifndef sphere_radius
DECLARE_PARAM(double, sphere_radius, 1.0)
#endif

//
// Boundary conditions
//
//- TODO: add description
#ifndef do_boundaries
DECLARE_PARAM(bool, do_boundaries, false)
#endif

//- TODO: add description
#ifndef stop_boundaries
DECLARE_PARAM(bool, stop_boundaries, false)
#endif

//- TODO: add description
#ifndef reflect_boundaries
DECLARE_PARAM(bool, reflect_boundaries, false)
#endif

#ifndef periodic_boundary_x
DECLARE_PARAM(bool, periodic_boundary_x, false)
#endif

#ifndef periodic_boundary_y
DECLARE_PARAM(bool, periodic_boundary_y, false)
#endif

#ifndef periodic_boundary_z
DECLARE_PARAM(bool, periodic_boundary_z, false)
#endif

//- tolerance to lattice mismatch for periodic boundaries:
//  when generating initial data with two lattice blocks
//  with different densities that need to be periodic,
//  allow this mismatch between lattice synchronization
//  at the boundary
//
#ifndef lattice_matchup_tolerance
DECLARE_PARAM(double, lattice_mismatch_tolerance, 0.05)
#endif

//
// I/O parameters
//
//- file prefix for input and intiial data file[s]
#ifndef initial_data_prefix
DECLARE_STRING_PARAM(initial_data_prefix, "initial_data")
#endif

//- ID-generator-specific parameter to overwrite initial data
#ifndef modify_initial_data
DECLARE_PARAM(bool, modify_initial_data, false)
#endif

//- file prefix for HDF5 output data file[s]
#ifndef output_h5data_prefix
DECLARE_STRING_PARAM(output_h5data_prefix, "output_data")
#endif

//- screen output frequency by iteration
#ifndef out_screen_every
DECLARE_PARAM(int32_t, out_screen_every, 1)
#endif

//- screen output frequency by time
//  if positive: overrides out_screen_every
#ifndef out_screen_dt
DECLARE_PARAM(double, out_screen_dt, 0)
#endif

//- scalar reductions output frequency by iteration
#ifndef out_scalar_every
DECLARE_PARAM(int32_t, out_scalar_every, 10)
#endif

//- scalar reductions output frequency by time
//  if positive: overrides out_scalar_every
#ifndef out_scalar_dt
DECLARE_PARAM(double, out_scalar_dt, 0)
#endif

// - diagnostic info output frequency by iteration
#ifndef out_diagnostic_every
DECLARE_PARAM(int32_t, out_diagnostic_every, 0);
#endif

//- HDF5 output frequency by iteration
#ifndef out_h5data_every
DECLARE_PARAM(int32_t, out_h5data_every, 10)
#endif

//- HDF5 output frequency by time
//  if positive: overrides out_h5data_every
#ifndef out_h5data_dt
DECLARE_PARAM(double, out_h5data_dt, 0)
#endif

//- produce separate HDF5 file per iteration
#ifndef out_h5data_separate_iterations
DECLARE_PARAM(bool, out_h5data_separate_iterations, false)
#endif

// WVT parameters
// Method:
// * Diehl et al., PASA 2015
// * Arth et al., 2019
//
// Boundary conditions:
// * reflective (default)
// * frozen
//
// wvt_mu
// Fraction of smoothing lenght that particles are
// allowed to move in one iteraiton
//
// wvt_ngb
// Number of desired particle neighbors
//
// wvt_convergence_check
// Stops simulation when convergence criterial is
// reached
//
// wvt_convergence_point
// Convergence criteria in terms of percentage of particles
// that moved more than 1.e-3 times the mean particle
// spacing at one iteration.
//
// wvt_h_ngb
// Alternative way to calculate the smoothing length based
// on the number of particle neighbors
//
// wvt_cool_down
// wvt iterations with decreasing wvt_mu
//
// wvt_radius
// sets wvt radius to apply boundary conditions

// - method for wvt pseudo-acceleration
#ifndef wvt_method
DECLARE_STRING_PARAM(wvt_method, "arth")
#endif

// - boundary condition for wvt particles
#ifndef wvt_boundary
DECLARE_STRING_PARAM(wvt_boundary, "reflective")
#endif

// - wvt_mu; control for step size
#ifndef wvt_mu
DECLARE_PARAM(double, wvt_mu, 1.e-3)
#endif

// - wvt_ngb; desired number of neighbors
#ifndef wvt_ngb
DECLARE_PARAM(double, wvt_ngb, 120)
#endif

// - switch for wvt convergence check
#ifndef wvt_convergence_check
DECLARE_PARAM(bool, wvt_convergence_check, true)
#endif

// - determined when simulation is converged
#ifndef wvt_convergence_point
DECLARE_PARAM(double, wvt_convergence_point, 1.0)
#endif

// - neighbor-based smoothing length
#ifndef wvt_h_ngb
DECLARE_PARAM(bool, wvt_h_ngb, false)
#endif

// - wvt iterations with decreasing step size
#ifndef wvt_cool_down
DECLARE_PARAM(int, wvt_cool_down, 0)
#endif

// - switch for wvt boundary setting
#ifndef wvt_set_boundary
DECLARE_PARAM(bool, wvt_set_boundary, true)
#endif

// - radius for wvt boundary conditions
#ifndef wvt_radius
DECLARE_PARAM(double, wvt_radius, 1.0)
#endif

//
// Viscosity and equation of state
//
//- which equation of state to use?
//  * "ideal fluid" (default)
//  * "polytropic"
//  * "white dwarf"
//  * "piecewise polytropic"
#ifndef eos_type
DECLARE_KEYWORD_PARAM(eos_type, eos_ideal)
#endif

// - file for tabulated EOS
#ifndef eos_tab_file_path
DECLARE_STRING_PARAM(eos_tab_file_path, ".")
#endif

//- polytropic index
#ifndef poly_gamma
DECLARE_PARAM(double, poly_gamma, 1.4)
#endif

//- additional polytropic index for piecewise polytrope
#ifndef poly_gamma2
DECLARE_PARAM(double, poly_gamma2, 2.5)
#endif

//- in piecewise polytropic equationa of state: threshold density
#ifndef ppt_density_thr
DECLARE_PARAM(double, ppt_density_thr, 5e+14)
#endif

// - defines viscosity prescription; options:
//   * constant: constant artificial_viscosity
//     cullen:   the Cullen'10 adaptive visc. prescription
#ifndef sph_viscosity
DECLARE_KEYWORD_PARAM(sph_viscosity,visc_constant)
#endif

//- for constant viscosity: parameter alpha (Rosswog'09, eq.59)
//  for adaptive Cullen+10 viscosity: minimum alpha
#ifndef sph_viscosity_alpha
DECLARE_PARAM(double, sph_viscosity_alpha, 1.0)
#endif

//- constant viscosity parameter beta (usually = 2*alpha)
#ifndef sph_viscosity_beta
DECLARE_PARAM(double, sph_viscosity_beta, 2.0)
#endif

//- viscosity parameter eta (avoids division by zero)
#ifndef sph_viscosity_epsilon
DECLARE_PARAM(double, sph_viscosity_epsilon, 0.01)
#endif

//- Cullen+10 viscosity: maximum allowed value of alpha
#ifndef sph_viscosity_alpha_max
  DECLARE_PARAM(double,sph_viscosity_alpha_max,1.0)
#endif

//- in adaptive Cullen+10 viscosity: how fast alpha returns to min alpha
//  tau = h / (2*l*v_sig)
//  typical values: 0.1 .. 0.2
#ifndef sph_viscosity_l
  DECLARE_PARAM(double,sph_viscosity_l,0.05)
#endif

//- in adaptive Cullen+10 viscosity: in the alpha_loc formula, relative 
//  weight between vsig^2 and A*h^2
#ifndef sph_viscosity_delta
  DECLARE_PARAM(double,sph_viscosity_delta,1.0)
#endif

//
// Gravity-related parameters
//
// Do FMM computation
#ifndef enable_fmm
DECLARE_PARAM(bool, enable_fmm, false)
#endif

//- mac'n'cheese acceptance criteria
#ifndef fmm_macangle
DECLARE_PARAM(double, fmm_macangle, 0.0)
#endif

//- maximum mass per cell
#ifndef fmm_max_cell_mass
DECLARE_PARAM(double, fmm_max_cell_mass, 0.)
#endif

//
// Parameters for particle relaxation, used to relax configurations
// by applying negative drag force against the direction of velocity
// for each particle:
//  f_relax = - (beta + gamma*v^2) * v
//
// Simple tests which are set up on regular rectangular lattices do not
// require particle relaxation term.
//

//- apply relaxation for this many steps (non-inclusive);
//  if set to zero (default), do not apply relaxation.
#ifndef relaxation_steps
DECLARE_PARAM(int, relaxation_steps, 0)
#endif

//- relaxation coefficients beta and gamma (both must be positive)
#ifndef relaxation_beta
DECLARE_PARAM(double, relaxation_beta, 1.e-6)
#endif

#ifndef relaxation_gamma
DECLARE_PARAM(double, relaxation_gamma, 0.0)
#endif

#ifndef relaxation_repulsion_radius
DECLARE_PARAM(double, relaxation_repulsion_radius, 0.25)
#endif

#ifndef relaxation_repulsion_gamma
DECLARE_PARAM(double, relaxation_repulsion_gamma, 0.0)
#endif

#ifndef evolve_internal_energy
DECLARE_PARAM(bool, evolve_internal_energy, true)
#endif

//
// Parameters for external acceleration
//
#ifndef thermokinetic_formulation
DECLARE_PARAM(bool, thermokinetic_formulation, true)
#endif

//- which external force to apply?
//  * "none" (default)
#ifndef external_force_type
DECLARE_STRING_PARAM(external_force_type, "none")
#endif

// poison zero potential level: since potential is defined
// up to a constant, any poison value should still work
#ifndef zero_potential_poison_value
DECLARE_PARAM(double, zero_potential_poison_value, 0.0)
#endif

// boundary wall power index
#ifndef extforce_wall_powerindex
DECLARE_PARAM(double, extforce_wall_powerindex, 5.0)
#endif

// boundary wall steepness parameter
#ifndef extforce_wall_steepness
DECLARE_PARAM(double, extforce_wall_steepness, 1e12)
#endif

// in mesa potential, fraction of the density-drop outer section to the radius
#ifndef mesa_rim_width
DECLARE_PARAM(double, mesa_rim_width, 0.25)
#endif

// gravitational acceleration constant on Earth
#ifndef gravity_acceleration_constant
DECLARE_PARAM(double, gravity_acceleration_constant, 9.81)
#endif

// Newtonian constant of gravitation (in CGS units)
#ifndef gravitational_constant
DECLARE_PARAM(double, gravitational_constant, 1)
//  DECLARE_PARAM(double,gravitational_constant, 6.674e-8)
#endif

//
// Parameters for the orbiting binary star setup
//
// binary orbital separation (in cm)
#ifndef orbital_separation
DECLARE_PARAM(double, orbital_separation, 2.5e9)
#endif

// in a NS-WD binary: mass of the neutron star (in g)
#ifndef mass_neutron_star
DECLARE_PARAM(double, mass_neutron_star, 1.26 * 1.988435e33)
#endif

// in a NS-WD binary: mass of the white dwarf (in g)
#ifndef mass_white_dwarf
DECLARE_PARAM(double, mass_white_dwarf, 1.10 * 1.988435e33)
#endif

//
// Specific apps
//
/// number of Sodtest to run (1..5)
#ifndef sodtest_num
DECLARE_PARAM(unsigned short, sodtest_num, 1)
#endif

// equal mass or equal particle separation switch for sodtube
#ifndef equal_mass
DECLARE_PARAM(bool, equal_mass, true)
#endif

// initial density profile: for initial data or density-supporting
// external portential
// Possible values:
// * 'constant'  :constant uniform-density spherical configuration
// * 'parabolic' :spherically-symmetric parabolic shape, rho ~ rho0*(1 - r^2)
// * 'mesa'      :constant density with a smooth parabolic fade-out on edge
// * 'from file' :setup density from the input_density_file
#ifndef density_profile
DECLARE_STRING_PARAM(density_profile, "constant")
#endif

// gridded input data for generating / supporting arbitrary density profiles
// used when parameter 'density_profile' is set to 'from file'
#ifndef input_density_file
DECLARE_STRING_PARAM(input_density_file, "")
#endif

// characteristic density for initial conditions
#ifndef rho_initial
DECLARE_PARAM(double, rho_initial, 1.0)
#endif

// characteristic pressure
#ifndef pressure_initial
DECLARE_PARAM(double, pressure_initial, 1.0)
#endif

// characteristic specific internal energy
#ifndef uint_initial
DECLARE_PARAM(double, uint_initial, 1.0)
#endif

// in Sedov test: total injected blast enregy
#ifndef sedov_blast_energy
DECLARE_PARAM(double, sedov_blast_energy, 1.0)
#endif

// in Sedov test: radius of energy injection
// (in units of particle separation)
#ifndef sedov_blast_radius
DECLARE_PARAM(double, sedov_blast_radius, 0.05)
#endif

// in Noh test: infall velocity
#ifndef noh_infall_velocity
DECLARE_PARAM(double, noh_infall_velocity, 0.1)
#endif

// initial data lattice type:
#ifndef lattice_type
DECLARE_PARAM(int, lattice_type, 0)
#endif

// if >0: lattice is randomly perturbed with this amplitude
// amplitude is in units of smoothing length (h)
#ifndef lattice_perturbation_amplitude
DECLARE_PARAM(double, lattice_perturbation_amplitude, 0.0)
#endif

// in several tests: initial velocity of the flow
#ifndef flow_velocity
DECLARE_PARAM(double, flow_velocity, 0.0)
#endif

// in several tests (e.g. KH and RT instabilities): density ratio
#ifndef density_ratio
DECLARE_PARAM(double, density_ratio, 2.0)
#endif

// A value from KH in Price's paper
#ifndef KH_A
DECLARE_PARAM(double, KH_A, 0.025)
#endif

// Lamdba value for KH in Price's paper
#ifndef KH_lambda
DECLARE_PARAM(double, KH_lambda, 1. / 6.)
#endif

// Rayleigh-Taylor instability: perturbation amplitude
#ifndef rt_perturbation_amplitude
DECLARE_PARAM(double, rt_perturbation_amplitude, 0.2)
#endif

// RT instability: the width of stripe where to apply perturbation
#ifndef rt_perturbation_stripe_width
DECLARE_PARAM(double, rt_perturbation_stripe_width, 0.1)
#endif

// RT instability: perturbation mode (1=one cusp, 2=two cusps etc.)
#ifndef rt_perturbation_mode
DECLARE_PARAM(double, rt_perturbation_mode, 1)
#endif

//
// Airfoil parameters
//
#ifndef airfoil_size
DECLARE_PARAM(double, airfoil_size, 2.0)
#endif

#ifndef airfoil_thickness
DECLARE_PARAM(double, airfoil_thickness, 0.05)
#endif

#ifndef airfoil_camber
DECLARE_PARAM(double, airfoil_camber, 0.1)
#endif

#ifndef airfoil_anchor_x
DECLARE_PARAM(double, airfoil_anchor_x, -1.0)
#endif

#ifndef airfoil_anchor_y
DECLARE_PARAM(double, airfoil_anchor_y, 0.0)
#endif

#ifndef airfoil_attack_angle
DECLARE_PARAM(double, airfoil_attack_angle, 0.0)
#endif

// ---

/*!
   trim whitespace from the begging and end of line
 */
std::string
trim(const std::string & str) {
  using namespace std;
  const size_t strBegin = str.find_first_not_of(" \t");
  if(strBegin == string::npos)
    return "";
  const size_t strEnd = str.find_last_not_of(" \t");
  return str.substr(strBegin, strEnd - strBegin + 1);
}

/**
 * @brief Sets a global parameter param_name to the value, given by its
 *        string representation param_value
 */
void
set_param(const std::string & param_name, const std::string & param_value) {

  // RANK/SIZE for LOG output
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  using namespace std;
  bool unknown_param = true;

  // strip trailing comments
  size_t i2 = param_value.find("#");
  string str_value;
  if(i2 != string::npos)
    str_value = trim(param_value.substr(0, i2));
  else
    str_value = param_value;

  // remove (single or double) quotes
  if(str_value[0] == '\'' or str_value[0] == '\"')
    str_value = str_value.substr(1, str_value.length() - 2);
  istringstream iss(str_value);

  // for boolean parameters
  bool lparam_value =
    (str_value == "yes" or str_value == "'yes'" or str_value == "\"yes\"" or
      str_value == "true" or str_value == "'true'" or str_value == "\"true\"");
  // timestepping and iterations --------------------------------------------
#ifndef initial_iteration
  READ_NUMERIC_PARAM(initial_iteration)
#endif

#ifndef final_iteration
  READ_NUMERIC_PARAM(final_iteration)
#endif

#ifndef initial_time
  READ_NUMERIC_PARAM(initial_time)
#endif

#ifndef final_time
  READ_NUMERIC_PARAM(final_time)
#endif

#ifndef initial_dt
  READ_NUMERIC_PARAM(initial_dt)
#endif

#ifndef timestep_cfl_factor
  READ_NUMERIC_PARAM(timestep_cfl_factor)
#endif

#ifndef adaptive_timestep
  READ_BOOLEAN_PARAM(adaptive_timestep)
#endif

# ifndef pressure_updates_number
  READ_NUMERIC_PARAM(pressure_updates_number)
# endif


  // particle number and density --------------------------------------------
#ifndef nparticles
  READ_NUMERIC_PARAM(nparticles)
#endif

#ifndef lattice_nx
  READ_NUMERIC_PARAM(lattice_nx)
#endif

#ifndef sph_eta
  READ_NUMERIC_PARAM(sph_eta)
#endif

#ifndef sph_smoothing_length
  READ_NUMERIC_PARAM(sph_smoothing_length)
#endif

#ifndef sph_separation
  READ_NUMERIC_PARAM(sph_separation)
#endif

  if(param_name == "sph_kernel") {
    for(int c = 0; c < str_value.length(); ++c)
      if(str_value[c] == ' ')
        str_value[c] = '_';

#ifndef sph_kernel
    if(boost::iequals(str_value, "cubic_spline"))
      _sph_kernel = cubic_spline;

    else if(boost::iequals(str_value, "quintic_spline"))
      _sph_kernel = quintic_spline;

    else if(boost::iequals(str_value, "wendland_c2"))
      _sph_kernel = wendland_c2;

    else if(boost::iequals(str_value, "wendland_c4"))
      _sph_kernel = wendland_c4;

    else if(boost::iequals(str_value, "wendland_c6"))
      _sph_kernel = wendland_c6;

    else if(boost::iequals(str_value, "gaussian"))
      _sph_kernel = gaussian;

    else if(boost::iequals(str_value, "super_gaussian"))
      _sph_kernel = super_gaussian;

    else if(boost::iequals(str_value, "sinc_ker"))
      _sph_kernel = sinc_ker;

    else {
      assert(false);
    }
#else
    if(not boost::iequals(str_value, QUOTE(sph_kernel))) {
      log_one(error) << "ERROR: sph_kernel #defined as \"" << QUOTE(sph_kernel)
                     << "\" "
                     << "but is reset to \"" << str_value
                     << "\" in parameter file" << std::endl;
      exit(2);
    }
#endif
    unknown_param = false;
  }

#ifndef sph_sinc_index
  READ_NUMERIC_PARAM(sph_sinc_index)
#endif

#ifndef sph_update_uniform_h
  READ_BOOLEAN_PARAM(sph_update_uniform_h)
#endif

#ifndef sph_variable_h
  READ_BOOLEAN_PARAM(sph_variable_h)
#endif

  // geometric configuration  -----------------------------------------------
#ifndef domain_type
  READ_NUMERIC_PARAM(domain_type)
#endif

#ifndef box_length
  READ_NUMERIC_PARAM(box_length)
#endif

#ifndef box_width
  READ_NUMERIC_PARAM(box_width)
#endif

#ifndef box_height
  READ_NUMERIC_PARAM(box_height)
#endif

#ifndef sphere_radius
  READ_NUMERIC_PARAM(sphere_radius)
#endif

  // boundary conditions  ---------------------------------------------------
#ifndef do_boundaries
  READ_BOOLEAN_PARAM(do_boundaries)
#endif

#ifndef stop_boundaries
  READ_BOOLEAN_PARAM(stop_boundaries)
#endif

#ifndef reflect_boundaries
  READ_BOOLEAN_PARAM(reflect_boundaries)
#endif

#ifndef periodic_boundary_x
  READ_BOOLEAN_PARAM(periodic_boundary_x)
#endif

#ifndef periodic_boundary_y
  READ_BOOLEAN_PARAM(periodic_boundary_y)
#endif

#ifndef periodic_boundary_z
  READ_BOOLEAN_PARAM(periodic_boundary_z)
#endif

#ifndef lattice_matchup_tolerance
  READ_NUMERIC_PARAM(lattice_mismatch_tolerance)
#endif

  // i/o parameters  --------------------------------------------------------
#ifndef initial_data_prefix
  READ_STRING_PARAM(initial_data_prefix)
#endif

#ifndef modify_initial_data
  READ_BOOLEAN_PARAM(modify_initial_data)
#endif

#ifndef output_h5data_prefix
  READ_STRING_PARAM(output_h5data_prefix)
#endif

#ifndef out_screen_every
  READ_NUMERIC_PARAM(out_screen_every)
#endif

#ifndef out_screen_dt
  READ_NUMERIC_PARAM(out_screen_dt)
#endif

#ifndef out_scalar_every
  READ_NUMERIC_PARAM(out_scalar_every)
#endif

#ifndef out_scalar_dt
  READ_NUMERIC_PARAM(out_scalar_dt)
#endif

#ifndef out_diagnostic_every
  READ_NUMERIC_PARAM(out_diagnostic_every)
#endif

#ifndef out_h5data_every
  READ_NUMERIC_PARAM(out_h5data_every)
#endif

#ifndef out_h5data_dt
  READ_NUMERIC_PARAM(out_h5data_dt)
#endif

#ifndef out_h5data_separate_iterations
  READ_BOOLEAN_PARAM(out_h5data_separate_iterations)
#endif

  // wvt parameters ---------------------------------------------------------
#ifndef wvt_method
  READ_STRING_PARAM(wvt_method)
#endif

#ifndef wvt_boundary
  READ_STRING_PARAM(wvt_boundary)
#endif

#ifndef wvt_mu
  READ_NUMERIC_PARAM(wvt_mu)
#endif

#ifndef wvt_ngb
  READ_NUMERIC_PARAM(wvt_ngb)
#endif

#ifndef wvt_convergence_check
  READ_BOOLEAN_PARAM(wvt_convergence_check)
#endif

#ifndef wvt_convergence_point
  READ_NUMERIC_PARAM(wvt_convergence_point)
#endif

#ifndef wvt_h_ngb
  READ_BOOLEAN_PARAM(wvt_h_ngb)
#endif

#ifndef wvt_cool_down
  READ_NUMERIC_PARAM(wvt_cool_down)
#endif

#ifndef wvt_radius
  READ_NUMERIC_PARAM(wvt_radius)
#endif

  // viscosity and equation of state ----------------------------------------
  if(param_name == "eos_type") {
    for(int c = 0; c < str_value.length(); ++c)
      if(str_value[c] == ' ')
        str_value[c] = '_';

#ifndef eos_type
    if(boost::iequals(str_value, "ideal_fluid"))
      _eos_type = eos_ideal;

    else if(boost::iequals(str_value, "polytropic"))
      _eos_type = eos_polytropic;

    else if(boost::iequals(str_value, "wd")
         or boost::iequals(str_value, "white_dwarf"))
      _eos_type = eos_wd;

    else if(boost::iequals(str_value, "ppt")
         or boost::iequals(str_value, "piecewise_polytropic"))
      _eos_type = eos_ppt;

    else if(boost::iequals(str_value, "no_eos")
         or boost::iequals(str_value, "none"))
      _eos_type = eos_no_eos;

    else {
      assert(false);
    }
#else
    if(not boost::iequals(str_value, QUOTE(eos_type))) {
      log_one(error) << "ERROR: eos_type #defined as \"" << QUOTE(eos_type)
                     << "\" "
                     << "but is reset to \"" << str_value
                     << "\" in parameter file" << std::endl;
      exit(2);
    }
#endif
    unknown_param = false;
  }

#ifndef eos_tab_file_path
  READ_STRING_PARAM(eos_tab_file_path)
#endif

#ifndef poly_gamma
  READ_NUMERIC_PARAM(poly_gamma)
#endif

#ifndef poly_gamma2
  READ_NUMERIC_PARAM(poly_gamma2)
#endif

#ifndef ppt_density_thr
  READ_NUMERIC_PARAM(ppt_density_thr)
#endif

// parsing sph_viscosity keywords
  if (param_name == "sph_viscosity") {
    for (int c=0; c<str_value.length(); ++c)
      if (str_value[c] == ' ') str_value[c] = '_';

#   ifndef sph_viscosity
    if (boost::iequals(str_value,"constant"))
      _sph_viscosity =       visc_constant;

    else if (boost::iequals(str_value,"cullen"))
      _sph_viscosity =            visc_cullen;

    else {
      log_one(error)
          << "ERROR: wrong value for sph_viscosity parameter"
          << std::endl;
      exit(2);
    }
#   else
    if (not boost::iequals(str_value,QUOTE(sph_viscosity))) {
      log_one(error)
          << "ERROR: sph_viscosity #define'd as \"" << QUOTE(sph_viscosity)
          << "\" but is reset to \"" << str_value << "\" in parameter file"
          << std::endl;
      exit(2);
    }
#   endif
    unknown_param = false;
  }

#ifndef sph_viscosity_alpha
  READ_NUMERIC_PARAM(sph_viscosity_alpha)
#endif

#ifndef sph_viscosity_beta
  READ_NUMERIC_PARAM(sph_viscosity_beta)
#endif

#ifndef sph_viscosity_epsilon
  READ_NUMERIC_PARAM(sph_viscosity_epsilon)
#endif

# ifndef sph_viscosity_alpha_max
  READ_NUMERIC_PARAM(sph_viscosity_alpha_max)
# endif

# ifndef sph_viscosity_l
  READ_NUMERIC_PARAM(sph_viscosity_l)
# endif

# ifndef sph_viscosity_delta
  READ_NUMERIC_PARAM(sph_viscosity_delta)
# endif

  // gravity-related  -------------------------------------------------------

#ifndef enable_fmm
  READ_BOOLEAN_PARAM(enable_fmm)
#endif

#ifndef fmm_macangle
  READ_NUMERIC_PARAM(fmm_macangle)
#endif

  // relaxation parameters  --------------------------------------------------
#ifndef relaxation_steps
  READ_NUMERIC_PARAM(relaxation_steps)
#endif

#ifndef relaxation_beta
  READ_NUMERIC_PARAM(relaxation_beta)
#endif

#ifndef relaxation_gamma
  READ_NUMERIC_PARAM(relaxation_gamma)
#endif

#ifndef relaxation_repulsion_radius
  READ_NUMERIC_PARAM(relaxation_repulsion_radius)
#endif

#ifndef relaxation_repulsion_gamma
  READ_NUMERIC_PARAM(relaxation_repulsion_gamma)
#endif

#ifndef evolve_internal_energy
  READ_BOOLEAN_PARAM(evolve_internal_energy)
#endif

  // external force  --------------------------------------------------------
#ifndef thermokinetic_formulation
  READ_BOOLEAN_PARAM(thermokinetic_formulation)
#endif

#ifndef external_force_type
  READ_STRING_PARAM(external_force_type)
#endif

#ifndef zero_potential_poison_value
  READ_NUMERIC_PARAM(zero_potential_poison_value)
#endif

#ifndef extforce_wall_powerindex
  READ_NUMERIC_PARAM(extforce_wall_powerindex)
#endif

#ifndef extforce_wall_steepness
  READ_NUMERIC_PARAM(extforce_wall_steepness)
#endif

#ifndef mesa_rim_width
  READ_NUMERIC_PARAM(mesa_rim_width)
#endif

#ifndef gravity_acceleration_constant
  READ_NUMERIC_PARAM(gravity_acceleration_constant)
#endif

#ifndef gravitational_constant
  READ_NUMERIC_PARAM(gravitational_constant)
#endif

#ifndef orbital_separation
  READ_NUMERIC_PARAM(orbital_separation)
#endif

#ifndef mass_neutron_star
  READ_NUMERIC_PARAM(mass_neutron_star)
#endif

#ifndef mass_white_dwarf
  READ_NUMERIC_PARAM(mass_white_dwarf)
#endif

  // specific apps  ---------------------------------------------------------
#ifndef sodtest_num
  READ_NUMERIC_PARAM(sodtest_num)
#endif

#ifndef equal_mass
  READ_BOOLEAN_PARAM(equal_mass)
#endif

#ifndef density_profile
  READ_STRING_PARAM(density_profile)
#endif

#ifndef input_density_file
  READ_STRING_PARAM(input_density_file)
#endif

#ifndef rho_initial
  READ_NUMERIC_PARAM(rho_initial)
#endif

#ifndef pressure_initial
  READ_NUMERIC_PARAM(pressure_initial)
#endif

#ifndef uint_initial
  READ_NUMERIC_PARAM(uint_initial)
#endif

#ifndef sedov_blast_energy
  READ_NUMERIC_PARAM(sedov_blast_energy)
#endif

#ifndef sedov_blast_radius
  READ_NUMERIC_PARAM(sedov_blast_radius)
#endif

#ifndef noh_infall_velocity
  READ_NUMERIC_PARAM(noh_infall_velocity)
#endif

#ifndef lattice_type
  READ_NUMERIC_PARAM(lattice_type)
#endif

#ifndef lattice_perturbation_amplitude
  READ_NUMERIC_PARAM(lattice_perturbation_amplitude)
#endif

#ifndef flow_velocity
  READ_NUMERIC_PARAM(flow_velocity)
#endif

#ifndef density_ratio
  READ_NUMERIC_PARAM(density_ratio)
#endif

#ifndef KH_A
  READ_NUMERIC_PARAM(KH_A)
#endif

#ifndef KH_lambda
  READ_NUMERIC_PARAM(KH_lambda)
#endif

#ifndef rt_perturbation_amplitude
  READ_NUMERIC_PARAM(rt_perturbation_amplitude)
#endif

#ifndef rt_perturbation_stripe_width
  READ_NUMERIC_PARAM(rt_perturbation_stripe_width)
#endif

#ifndef rt_perturbation_mode
  READ_NUMERIC_PARAM(rt_perturbation_mode)
#endif

  // airfoil parameters  ----------------------------------------------------
#ifndef airfoil_size
  READ_NUMERIC_PARAM(airfoil_size)
#endif

#ifndef airfoil_thickness
  READ_NUMERIC_PARAM(airfoil_thickness)
#endif

#ifndef airfoil_camber
  READ_NUMERIC_PARAM(airfoil_camber)
#endif

#ifndef airfoil_anchor_x
  READ_NUMERIC_PARAM(airfoil_anchor_x)
#endif

#ifndef airfoil_anchor_y
  READ_NUMERIC_PARAM(airfoil_anchor_y)
#endif

#ifndef airfoil_attack_angle
  READ_NUMERIC_PARAM(airfoil_attack_angle)
#endif

  // unknown parameter -------------------------------
  if(unknown_param) {
    log_one(error) << "ERROR: unknown parameter " << param_name << endl;
    exit(2);
  }

  log_one(trace) << param_name << ": " << param_value << endl;
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
void
read_params(const char * parfile) {
  using namespace std;
  ifstream infile;
  string line;

  // attempt to open the file
  infile.open(parfile);
  if(!infile) {
    cerr << "ERROR: Unable to open parameter file 'input.par'" << endl;
    exit(1);
  }

  for(int ln = 1; std::getline(infile, line); ++ln) {

    // skip comments (lines starting with '#' at any position)
    // and blank lines
    bool is_blank = true, is_comment = false;
    for(size_t i = 0; i < line.length(); i++) {
      char c = line[i];
      is_comment = (c == '#');
      is_blank = (c == ' ' || c == '\t');
      if(is_comment or not is_blank)
        break;
    }
    if(is_comment or is_blank)
      continue;

    // line must be in the form:
    // lhs = rhs
    size_t eqpos = line.find("=");
    if(eqpos == string::npos) {
      cerr << "ERROR in parameter file " << parfile << ", line #" << ln << endl;
      exit(1);
    }
    string lhs = line.substr(0, eqpos);
    string rhs = line.substr(eqpos + 1);

    // strip whitespace
    lhs = trim(lhs);
    rhs = trim(rhs);
    if(lhs.length() == 0 or rhs.length() == 0) {
      cerr << "ERROR in parameter file " << parfile << ", line #" << ln
           << ": wrong format" << endl;
      exit(1);
    }
    set_param(lhs, rhs);
  }
  infile.close();
}

/**
 * @brief MPI parameter file reader
 * @todo  Use FleCSI infrastructure instead (e.g. Flecsi_Sim_IO?)
 */
void
mpi_read_params(const char * parameter_file) {
  const int MAXLEN = 2048;
  char buffer[MAXLEN];
  char * parfile = buffer;
  int len, rank, size, parfile_free = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;

  // the char* array pointer is only valid on rank 0;
  // broadcast parfile name from rank 0 to other ranks
  if(rank == 0) {
    strcpy(parfile, parameter_file);
    len = strlen(parfile);
  }
  MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(parfile, len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

  log_one(trace) << "Parameter file name on rank " << rank << " over " << size
                 << ": " << parfile << std::endl
                 << std::flush;

  // queue ranks to read the parfile sequentially;
  // wait for a message from previous rank, unless this is rank 0
  if(rank > 0)
    MPI_Recv(&parfile_free, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &status);

  // read parameters
  read_params(parfile);

  // inform the next rank that the file is free to read
  if(rank + 1 < size)
    MPI_Send(&parfile_free, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);

} // mpi_read_param

} // namespace param

#endif // PARAMS_H
