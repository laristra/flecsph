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

#include "analysis.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "diagnostic.h"

#define OUTPUT_ANALYSIS

static std::string output_h5data_file; // = output_h5data_prefix + ".h5part"

using namespace flecsph_log;

void
set_derived_params() {
  using namespace param;

  // set kernel
  kernels::select();

  // set viscosity
  viscosity::select();

  // filenames (this will change for multiple files output)
  std::ostringstream oss;
  oss << output_h5data_prefix << ".h5part";
  output_h5data_file = oss.str();

  // analysis: set output times
  analysis::set_initial_time_iteration();

  // set equation of state
  eos::select();

  // set gravitational constant
  fmm::gc = gravitational_constant;

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

  // set simulation parameters
  param::mpi_read_params(parameter_file);
  set_derived_params();

  // read input file and initialize equation of state
  body_system<double, gdimension> bs;
  bs.read_bodies(initial_data_prefix, output_h5data_prefix, initial_iteration);
  bs.setMacangle(param::fmm_macangle);

  MPI_Barrier(MPI_COMM_WORLD);

  do {
    analysis::screen_output(rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if(physics::iteration == param::initial_iteration) {

      log_one(trace) << "First iteration" << std::endl;
      bs.update_iteration();
      bs.apply_all(eos::init);

      if (sph_viscosity != visc_constant) {
        bs.apply_all(viscosity::initialize_alpha);
      }

      if(thermokinetic_formulation) {
        // at this point, gravitational potential is not set yet
        // we set total energy nevertheless, because thermodynamic
        // quantities (rho, P, cs) need to be computed still
        bs.apply_all(physics::set_total_energy);
      }

      log_one(trace) << "compute density pressure cs" << std::endl;
      bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
      bs.apply_all(integration::save_velocityhalf);

      if (sph_viscosity != visc_constant) {
        log_one(trace) << "computing adaptive viscosity" << std::endl;
        bs.apply_in_smoothinglength(viscosity::compute_alpha);
      }

      // compute acceleration
      log_one(trace) << "compute rhs of evolution equations" << std::endl;
      bs.reset_ghosts();
      bs.apply_in_smoothinglength(physics::compute_acceleration);
      if(param::enable_fmm){
        log_one(trace) << "compute gravitation" << std::endl;
        bs.gravitation_fmm();
      }
      if (physics::iteration < relaxation_steps) {
        log_one(trace) << "add relaxation terms" << std::endl;
        bs.apply_all(physics::add_drag_acceleration);
        bs.apply_in_smoothinglength(physics::add_short_range_repulsion);
        log_one(trace) << "relaxation terms: done" << std::endl;
      }

      if(thermokinetic_formulation) {
        // compute total energy for every particle
        bs.apply_all(physics::set_total_energy);
      }

      if (adaptive_timestep) {
        // Update timestep in the very beginning
        log_one(trace) << "compute adaptive timestep" << std::endl;
        bs.apply_all(physics::compute_dt);
        bs.get_all(physics::set_adaptive_timestep);
        log_one(trace) << ".done" << std::endl;
      }

      if (evolve_internal_energy) {
        if (thermokinetic_formulation){
          // compute de/dt 
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dedt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dedt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dedt);

            bs.apply_all(physics::recompute_pressure_soundspeed_thermokinetic);
            if (m < pressure_updates_number) 
              bs.reset_ghosts(); // skip syncing with the last pass
          }
        }else{
          // or compute du/dt
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dudt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dudt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dudt);
            bs.apply_all(physics::recompute_pressure_soundspeed);
            if (m < pressure_updates_number) 
              bs.reset_ghosts(); // skip syncing with the last pass
          }
        }
      } // if evolve_internal_energy
      log_one(trace) << "compute initial rhs terms: done" << std::endl;
    }
    else { // not the initial iteration
      log_one(trace) << "leapfrog: kick one" << std::endl;
      if (evolve_internal_energy) {
        if (thermokinetic_formulation)
          bs.apply_all(integration::leapfrog_kick_e);
        else
          bs.apply_all(integration::leapfrog_kick_u);
      }
      bs.apply_all(integration::leapfrog_kick_v);
      bs.apply_all(integration::save_velocityhalf);
      log_one(trace) << "kick one: done" << std::endl;

      log_one(trace) << "leapfrog: drift" << std::endl;
      bs.apply_all(integration::leapfrog_drift);
      log_one(trace) << "drift: done" << std::endl;

      // sync velocities
      bs.update_iteration();
      log_one(trace) << "compute density pressure cs" << std::endl;
      bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);

      if (sph_viscosity != visc_constant) {
        log_one(trace) << "computing adaptive viscosity" << std::endl;
        bs.apply_in_smoothinglength(viscosity::compute_alpha);
      }

      // compute acceleration
      log_one(trace) << "leapfrog: kick two (velocity)" << std::endl;
      bs.reset_ghosts();
      bs.apply_in_smoothinglength(physics::compute_acceleration);
      if(param::enable_fmm){
        log_one(trace) << "computing gravitation" << std::endl;
        bs.gravitation_fmm();
      }
      if (physics::iteration < relaxation_steps) {
        bs.apply_all(physics::add_drag_acceleration);
        bs.apply_in_smoothinglength(physics::add_short_range_repulsion);
      }
      bs.apply_all(integration::leapfrog_kick_v);
      log_one(trace) << "kick two (velocity): done" << std::endl;

      // sync velocities: needed for de/dt (du/dt)
      bs.reset_ghosts();

      if (evolve_internal_energy) {
        log_one(trace) << "leapfrog: kick two (energy)" << std::endl;
        if (thermokinetic_formulation) {
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dedt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dedt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dedt);

            bs.apply_all(physics::recompute_pressure_soundspeed_thermokinetic);
            if (m < pressure_updates_number)
              bs.reset_ghosts(); // skip syncing with the last pass
          }

          bs.apply_all(integration::leapfrog_kick_e);
        }
        else {
          for (int m=1; m<=pressure_updates_number;++m) { // 1 or 2 passes
            log_one(trace) << "compute dudt: pass " << m  << std::endl;
            bs.apply_in_smoothinglength(physics::compute_dudt);
            if (physics::iteration < relaxation_steps)
              bs.apply_all(physics::add_drag_dudt);

            bs.apply_all(physics::recompute_pressure_soundspeed);
            if (m < pressure_updates_number) 
              bs.reset_ghosts(); // skip syncing with the last pass
          }
          bs.apply_all(integration::leapfrog_kick_u);
        }
        log_one(trace) << "kick two (energy): done" << std::endl;
      } // evolve internal energy
    } // not initial iteration

    if(sph_variable_h){
      log_one(trace) << "updating smoothing length"<<std::endl;
      bs.get_all(physics::compute_smoothinglength);
      log_one(trace) << ".done" << std::endl;
    }else if(sph_update_uniform_h){
      // The particles moved, compute new smoothing length
      log_one(trace) << "updating smoothing length"<<std::endl;
      bs.get_all(physics::compute_average_smoothinglength,bs.getNBodies());
      log_one(trace) << ".done" << std::endl;
    }

    // Compute and output scalar reductions and diagnostic
    analysis::scalar_output(bs,rank);
    analysis::h5data_output(bs, rank);
    diagnostic::output(bs,rank);

    // Check for nans
    bs.apply_all(physics::check_nans);
    bs.apply_all(physics::check_negativity);

    if(adaptive_timestep) {
      // Update timestep
      log_one(trace) << "compute adaptive timestep" << std::endl;
      bs.apply_all(physics::compute_dt);
      bs.get_all(physics::set_adaptive_timestep);
      log_one(trace) << ".done" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    physics::advance_time();

  } while(not physics::termination_criteria());
} // mpi_init_task

flecsi_register_mpi_task(mpi_init_task, flecsi::execution);

void
usage() {
  log_one(warn) << "Usage: ./hydro_" << gdimension << "d "
                    << "<parameter-file.par>" << std::endl;
}

bool
check_conservation(const std::vector<analysis::e_conservation> & check) {
  return analysis::check_conservation(check);
}

void
specialization_tlt_init(int argc, char * argv[]) {
  log_set_output_rank(0);

  log_one(trace) << "In user specialization_driver" << std::endl;

  // check options list: exactly one option is allowed
  if(argc != 2) {
    log_one(error) << "ERROR: parameter file not specified!" << std::endl;
    usage();
    return;
  }

  flecsi_execute_mpi_task(mpi_init_task, flecsi::execution, argv[1]);

} // specialization driver

void
driver(int, char **) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  log_one(trace) << "In user driver" << std::endl;
} // driver

} // namespace execution
} // namespace flecsi
