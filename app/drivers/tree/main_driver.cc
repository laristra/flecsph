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

  // filenames (this will change for multiple files output)
  std::ostringstream oss;
  oss << output_h5data_prefix << ".h5part";
  output_h5data_file = oss.str();

  // iteration and time
  physics::iteration = initial_iteration;
  physics::totaltime = initial_time;
  physics::dt = initial_dt;

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
      
  size_t total = 5; 
  do {
    analysis::screen_output(rank);
    bs.update_iteration();
    double begin = omp_get_wtime(); 
    size_t total = 0; 
    bs.apply_in_smoothinglength(
      [&](tree_topology_t::entity_t& e, std::vector<tree_topology_t::entity_t*> & n, size_t& total){
        total+=n.size(); 
      },total
    );
    std::cout<<"Average: "<<total/bs.getNBodies()<<std::endl;
    double end = omp_get_wtime(); 
    std::cout<<"Traversal time: "<<end-begin<<"s "<<std::endl;
    ++physics::iteration;
  } while(--total!=0);
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