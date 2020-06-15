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
 * @file main_driver.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Specialization and Main driver used in FleCSI.
 * The Specialization Driver is normally used to register data and the main
 * code is in the Driver.
 */

#include <iostream>
#include <numeric> // For accumulate

#include <mpi.h>
#ifdef ENABLE_LEGION
#include <legion.h>
#endif
#include <omp.h>

#include "flecsi/data/data.h"
#include "flecsi/data/data_client.h"
#include "flecsi/execution/execution.h"

// #define poly_gamma 5./3.
#include "analysis.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "diagnostic.h"
#include "params.h"
#include "wvt.h"

#define OUTPUT_ANALYSIS

static std::string output_h5data_file; // = output_h5data_prefix + ".h5part"

void
set_derived_params() {
  using namespace param;

  // set kernel
  kernels::select();

  // set viscosity
  viscosity::select(sph_viscosity);

  // density profiles
  density_profiles::select();

  // wvt method
  wvt::select();

  // filenames (this will change for multiple files output)
  std::ostringstream oss;
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

namespace flecsi {
namespace execution {

void
mpi_init_task(const char * parameter_file) {
  using namespace param;

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  log_one(info) << "" << std::endl;
  log_one(info) << "  ---------------------------------------------- "
                << std::endl;
  log_one(info) << " |          Multi-D WVT Relaxation Driver       |"
                << std::endl;
  log_one(info) << "  ---------------------------------------------- "
                << std::endl;
  log_one(info) << "! Caution: The driver is currently only working !"
                << std::endl;
  log_one(info) << "!   for objects that are SPHERICALLY SYMMETRIC  !"
                << std::endl;
  log_one(info) << "" << std::endl;
  log_one(info) << "" << std::endl;

  // set simulation parameters
  param::mpi_read_params(parameter_file);
  set_derived_params();

  // read input file and initialize equation of state
  body_system<double, gdimension> bs;
  bs.read_bodies(initial_data_prefix, output_h5data_prefix, initial_iteration);

  MPI_Barrier(MPI_COMM_WORLD);

  do {
    analysis::screen_output(rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if(physics::iteration == param::initial_iteration) {
      log_one(trace) << "First iteration" << std::endl << std::flush;
      bs.update_iteration();
      bs.apply_all(eos::init);

      log_one(trace) << "compute density (for output)" << std::endl
                     << std::flush;
      bs.apply_in_smoothinglength(wvt::compute_density);

      // necessary?
      bs.reset_ghosts();

      log_one(trace) << "compute wvt acceleration" << std::endl << std::flush;
      bs.apply_in_smoothinglength(wvt::wvt_acceleration);
      log_one(trace) << ".done" << std::endl;
    }
    else if(physics::iteration <= final_iteration) {
      log_one(trace) << "wvt displacement" << std::flush;
      bs.apply_all(wvt::wvt_displacement);
      log_one(trace) << ".done" << std::endl;

      if(sph_variable_h) {
        log_one(trace) << "updating wvt smoothing length" << std::flush;
        bs.get_all(wvt::compute_smoothinglength_wvt);
        log_one(trace) << ".done" << std::endl << std::flush;
      }

      // sync velocities?
      bs.update_iteration();
      log_one(trace) << "compute density (for output)" << std::endl
                     << std::flush;
      bs.apply_in_smoothinglength(wvt::compute_density);
      log_one(trace) << ".done" << std::endl;

      bs.get_all(wvt::calculate_standard_deviation);

      // necessary?
      bs.reset_ghosts();

      log_one(trace) << "compute wvt acceleration" << std::endl << std::flush;
      bs.apply_in_smoothinglength(wvt::wvt_acceleration);
      log_one(trace) << ".done" << std::endl;

      // sync velocities
      bs.reset_ghosts();

      if((physics::iteration > param::initial_iteration) &&
         wvt_convergence_check) {
        log_one(trace) << "check convergence" << std::endl << std::flush;
        bs.get_all(wvt::check_convergence_wvt);
        log_one(trace) << ".done" << std::endl;
      }
    }
    else {
      log_one(trace) << "wvt cool down" << std::endl;

      log_one(trace) << "wvt displacement" << std::flush;
      bs.apply_all(wvt::wvt_displacement);
      log_one(trace) << ".done" << std::endl;

      if(sph_variable_h) {
        log_one(trace) << "updating wvt smoothing length" << std::flush;
        bs.get_all(wvt::compute_smoothinglength_wvt);
        log_one(trace) << ".done" << std::endl << std::flush;
      }

      bs.update_iteration();
      log_one(trace) << "compute density (for output)" << std::endl
                     << std::flush;
      bs.apply_in_smoothinglength(wvt::compute_density);
      log_one(trace) << ".done" << std::endl;

      bs.get_all(wvt::calculate_standard_deviation);

      bs.reset_ghosts();

      log_one(trace) << "compute wvt acceleration" << std::endl << std::flush;
      bs.apply_in_smoothinglength(wvt::wvt_acceleration);
      log_one(trace) << ".done" << std::endl;

      bs.reset_ghosts();
    }

    // Compute and output scalar reductions and diagnostic
    analysis::scalar_output(bs, rank);
    diagnostic::output(bs, rank);

    if((wvt_basic::wvt_converged) ||
       (physics::iteration == final_iteration + wvt_cool_down)) {
      log_one(trace) << "reset density to profile" << std::endl << std::flush;
      bs.apply_all(wvt::wvt_set_density);
      log_one(trace) << ".done" << std::endl;
    }

    if((out_h5data_every > 0 && physics::iteration % out_h5data_every == 0) ||
       (wvt_basic::wvt_converged == true) ||
       (physics::iteration == (final_iteration + wvt_cool_down))) {
      bs.write_bodies(
        output_h5data_prefix, physics::iteration, physics::totaltime);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    ++physics::iteration;

    physics::totaltime += physics::dt;

  } while(physics::iteration <= (final_iteration + wvt_cool_down) &&
          wvt_basic::wvt_converged == false);
} // mpi_init_task

flecsi_register_mpi_task(mpi_init_task, flecsi::execution);

void
usage(int rank) {
  log_one(warn) << "Usage: ./hydro_" << gdimension << "d "
                << "<parameter-file.par>" << std::endl
                << std::flush;
}

bool

check_conservation(const std::vector<analysis::e_conservation> & check) {
  return analysis::check_conservation(check);
}

void
specialization_tlt_init(int argc, char * argv[]) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  log_set_output_rank(0);

  log_one(trace) << "In user specialization_driver" << std::endl;

  // check options list: exactly one option is allowed
  if(argc != 2) {
    log_one(error) << "ERROR: parameter file not specified!" << std::endl;
    usage(rank);
    return;
  }

  flecsi_execute_mpi_task(mpi_init_task, flecsi::execution, argv[1]);

} // specialization driver

void
driver(int argc, char * argv[]) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  log_one(trace) << "In user driver" << std::endl;
} // driver

} // namespace execution
} // namespace flecsi
