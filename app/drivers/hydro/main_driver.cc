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
 * @file main_driver.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Specialization and Main driver used in FleCSI.
 * The Specialization Driver is normally used to register data and the main
 * code is in the Driver.
 */

#include <iostream>
#include <numeric> // For accumulate
#include <iostream>

#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

// #define poly_gamma 5./3.
#include "params.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "analysis.h"

#define OUTPUT_ANALYSIS

static std::string initial_data_file;  // = initial_data_prefix  + ".h5part"
static std::string output_h5data_file; // = output_h5data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  // set kernel
  kernels::select(sph_kernel);

  // filenames (this will change for multiple files output)
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();
  oss << output_h5data_prefix << ".h5part";
  output_h5data_file = oss.str();

  // iteration and time
  physics::iteration = initial_iteration;
  physics::totaltime = initial_time;
  physics::dt = initial_dt; // TODO: use particle separation and Courant factor

  // set equation of state
  eos::select(eos_type);

  // set external force
  external_force::select(external_force_type);
}

namespace flecsi{
namespace execution{

void
mpi_init_task(const char * parameter_file){
  using namespace param;

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // set simulation parameters
  param::mpi_read_params(parameter_file);
  set_derived_params();

  // remove output file
  remove(output_h5data_file.c_str());

  // read input file
  body_system<double,gdimension> bs;
  bs.read_bodies(initial_data_file.c_str(),initial_iteration);

  // boundaries
/*
  auto range_boundaries = bs.getRange();
  point_t distance = range_boundaries[1]-range_boundaries[0];
  for(unsigned short i = 0; i < gdimension; ++i){
    distance[i] = fabs(distance[i]);
  }
  double h = bs.getSmoothinglength();
  physics::min_boundary = {(0.1+2*h)*distance+range_boundaries[0]};
  physics::max_boundary = {-(0.1-2*h)*distance+range_boundaries[1]};
  clog_one(info) << "Limits: " << physics::min_boundary << " ; "
         << physics::max_boundary << std::endl;
 */
  MPI_Barrier(MPI_COMM_WORLD);
  bs.update_iteration();
  if(thermokinetic_formulation) {
    // compute total energy for every particle
    bs.apply_all(physics::set_total_energy);
  }

  bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
    
  if(out_scalar_every > 0 && physics::iteration % out_scalar_every == 0){
    // Compute conserved quantities
    bs.get_all(analysis::compute_lin_momentum);
    bs.get_all(analysis::compute_total_mass);
    bs.get_all(analysis::compute_total_energy);
    bs.get_all(analysis::compute_total_ang_mom);
    analysis::scalar_output("scalar_reductions.dat");
  }

  bs.write_bodies(output_h5data_prefix,physics::iteration);

  ++physics::iteration;
  do {
    analysis::screen_output(rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if (physics::iteration == 1){
      // Compute and prepare the tree for this iteration
      // - Compute the Max smoothing length
      // - Compute the range of the system using the smoothinglength
      // - Compute the keys
      // - Distributed qsort and sharing
      // - Generate and feed the tree
      // - Exchange branches for smoothing length
      // - Compute and exchange ghosts in real smoothing length
      bs.update_iteration();

      // at the initial iteration, P, rho and cs have not been computed yet;
      // for all subsequent steps, however, they are computed at the end 
      // of the iteration
      rank|| clog(trace) << "first iteration: pressure, rho and cs" << std::flush;
      bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
      bs.apply_all(physics::save_velocityhalf);
      rank|| clog(trace) << ".done" << std::endl;

      // necessary for computing dv/dt and du/dt in the next step
      bs.update_neighbors();

      rank|| clog(trace) << "compute rhs of evolution equations" << std::flush;
      bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
      if (thermokinetic_formulation)
        bs.apply_in_smoothinglength(physics::compute_dedt);
      else
        bs.apply_in_smoothinglength(physics::compute_dudt);
      rank|| clog(trace) << ".done" << std::endl;

    }
    else {
      rank|| clog(trace) << "leapfrog: kick one" << std::flush;
      bs.apply_all(physics::leapfrog_kick_v);
      if (thermokinetic_formulation)
        bs.apply_all(physics::leapfrog_kick_e);
      else
        bs.apply_all(physics::leapfrog_kick_u);
      bs.apply_all(physics::save_velocityhalf);
      rank|| clog(trace) << ".done" << std::endl;

      // sync velocities
      bs.update_neighbors();

      rank|| clog(trace) << "leapfrog: drift" << std::flush;
      bs.apply_all(physics::leapfrog_drift);
      bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
      rank|| clog(trace) << ".done" << std::endl;

      // recompute the tree and sync
      bs.update_iteration();
      bs.update_neighbors();

      rank|| clog(trace) << "leapfrog: kick two (velocity)" << std::flush;
      bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
      bs.apply_all(physics::leapfrog_kick_v);
      rank|| clog(trace) << ".done" << std::endl;

      // sync velocities
      bs.update_neighbors();

      rank|| clog(trace) << "leapfrog: kick two (energy)" << std::flush;
      if (thermokinetic_formulation) {
        bs.apply_in_smoothinglength(physics::compute_dedt);
        bs.apply_all(physics::leapfrog_kick_e);
      }
      else {
        bs.apply_in_smoothinglength(physics::compute_dudt);
        bs.apply_all(physics::leapfrog_kick_u);
      }
      rank|| clog(trace) << ".done" << std::endl;
    }

    if(sph_update_uniform_h){
      // The particles moved, compute new smoothing length 
      rank || clog(trace) << "updating smoothing length"<<std::flush;
      bs.get_all(physics::compute_average_smoothinglength,bs.getNBodies());
      rank || clog(trace) << ".done" << std::endl << std::flush;
    }


    if(out_scalar_every > 0 && physics::iteration % out_scalar_every == 0){
      // Compute the analysis values based on physics
      bs.get_all(analysis::compute_lin_momentum);
      bs.get_all(analysis::compute_total_mass);
      bs.get_all(analysis::compute_total_energy);
      bs.get_all(analysis::compute_total_ang_mom);
      // Only add the header in the first iteration
      analysis::scalar_output("scalar_reductions.dat");
    }


#ifdef OUTPUT
    if(out_h5data_every > 0 && physics::iteration % out_h5data_every == 0){
      bs.write_bodies(output_h5data_prefix,physics::iteration/out_h5data_every);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    ++physics::iteration;
    physics::totaltime += physics::dt;

  } while(physics::iteration <= final_iteration);
} // mpi_init_task


flecsi_register_mpi_task(mpi_init_task, flecsi::execution);

void 
usage(int rank) {
  rank|| clog(warn) << "Usage: ./hydro_" << gdimension << "d " 
                    << "<parameter-file.par>" << std::endl << std::flush;
}


void
specialization_tlt_init(int argc, char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  rank|| clog(trace) << "In user specialization_driver" << std::endl;

  // check options list: exactly one option is allowed
  if (argc != 2) {
    rank|| clog(error) << "ERROR: parameter file not specified!" << std::endl;
    usage(rank);
    return;
  }

  flecsi_execute_mpi_task(mpi_init_task, flecsi::execution, argv[1]);

} // specialization driver


void
driver(int argc,  char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  rank|| clog(trace) << "In user driver" << std::endl;
} // driver


} // namespace execution
} // namespace flecsi


