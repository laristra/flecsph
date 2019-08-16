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
#include <iostream>

#include <mpi.h>
#ifdef ENABLE_LEGION
#include <legion.h>
#endif
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

// #define poly_gamma 5./3.
#include "params.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "wvt.h"
#include "analysis.h"
#include "diagnostic.h"

#define OUTPUT_ANALYSIS

static std::string output_h5data_file; // = output_h5data_prefix + ".h5part"

void set_derived_params() {
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




  // read input file and initialize equation of state
  body_system<double,gdimension> bs;
  bs.read_bodies(initial_data_prefix,
      output_h5data_prefix,initial_iteration);

  MPI_Barrier(MPI_COMM_WORLD);

  do {
    analysis::screen_output(rank);
    MPI_Barrier(MPI_COMM_WORLD);

    if (physics::iteration == param::initial_iteration){

      clog_one(trace)<<"First iteration"<<std::endl << std::flush;
      bs.update_iteration();
      bs.apply_all(eos::init);

      clog_one(trace) << "compute density (for output)"<<std::endl << std::flush;
      bs.apply_in_smoothinglength(physics::compute_density);

      // necessary?
      bs.reset_ghosts();

      clog_one(trace) << "compute wvt acceleration"<<std::endl << std::flush;
      bs.apply_in_smoothinglength(wvt::wvt_acceleration);
      clog_one(trace) << ".done" << std::endl;
    }
    else {
      clog_one(trace) << "wvt displacement" << std::flush;
      bs.apply_all(wvt::wvt_displacement);
      clog_one(trace) << ".done" << std::endl;

      // sync velocities?
      bs.update_iteration();
      clog_one(trace) << "compute density (for output)"<<std::endl << std::flush;
      bs.apply_in_smoothinglength(physics::compute_density);
      clog_one(trace) << ".done" << std::endl;

      // necessary?
      bs.reset_ghosts();

      clog_one(trace) << "compute wvt acceleration"<<std::endl << std::flush;
      bs.apply_in_smoothinglength(wvt::wvt_acceleration);
      clog_one(trace) << ".done" << std::endl;

      // sync velocities
      bs.reset_ghosts();
    }

    if(sph_variable_h){
     clog_one(trace) << "updating wvt smoothing length"<<std::flush;
      bs.get_all(wvt::compute_smoothinglength_wvt);
      clog_one(trace) << ".done" << std::endl << std::flush;
    }

//    if (adaptive_timestep) {
//      // Update timestep
//      clog_one(trace) << "compute adaptive timestep" << std::flush;
//      bs.apply_in_smoothinglength(physics::estimate_maxmachnumber);
//      bs.apply_all(physics::compute_dt);
//      bs.get_all(physics::set_adaptive_timestep);
//      clog_one(trace) << ".done" << std::endl;
//    }

    // Compute and output scalar reductions and diagnostic
    analysis::scalar_output(bs,rank);
    diagnostic::output(bs,rank);

    if(out_h5data_every > 0 && physics::iteration % out_h5data_every == 0){
      bs.write_bodies(output_h5data_prefix,physics::iteration,
          physics::totaltime);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ++physics::iteration;

    physics::totaltime += physics::dt;

  } while(physics::iteration <= final_iteration);
} // mpi_init_task


flecsi_register_mpi_task(mpi_init_task, flecsi::execution);

void
usage(int rank) {
  clog_one(warn) << "Usage: ./hydro_" << gdimension << "d "
                    << "<parameter-file.par>" << std::endl << std::flush;
}

bool
check_conservation(
  const std::vector<analysis::e_conservation>& check
)
{
  return analysis::check_conservation(check);
}

void
specialization_tlt_init(int argc, char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  clog_set_output_rank(0);

  clog_one(trace) << "In user specialization_driver" << std::endl;

  // check options list: exactly one option is allowed
  if (argc != 2) {
    clog_one(error) << "ERROR: parameter file not specified!" << std::endl;
    usage(rank);
    return;
  }

  flecsi_execute_mpi_task(mpi_init_task, flecsi::execution, argv[1]);

} // specialization driver


void
driver(int argc,  char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  clog_one(trace) << "In user driver" << std::endl;
} // driver


} // namespace execution
} // namespace flecsi
