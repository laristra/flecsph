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

#include "bodies_system.h"
#include "default_physics.h"
#include "analysis.h"

#define OUTPUT_ANALYSIS

namespace flecsi{
namespace execution{

void
mpi_init_task(int numberiterations){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  clog_set_output_rank(0);
  
  int iteroutput = 1;
  double maxtime = 10.0; // TODO: this is never used.

  // Init if default values are not ok
  physics::dt = 0.0025;
  physics::alpha = 1; 
  physics::beta = 2; 
  physics::do_boundaries = true;
  physics::stop_boundaries = true;
  physics::gamma = 1.4;
  physics::epsilon = 0.01;

  const char * inputFile = "hdf5_sodtube.h5part";
  const char * outputFile = "output_sodtube.h5part"; 
  // Remove output file 
  remove(outputFile); 

  body_system<double,gdimension> bs;
  bs.read_bodies("hdf5_sodtube.h5part",physics::iteration);

  // Set the boundaries to be at 10% of the total range
  double h = bs.getSmoothinglength();
  auto range_boundaries = bs.getRange(); 
  double distance = fabs(range_boundaries[1][0]-range_boundaries[0][0]);
  physics::min_boundary = {range_boundaries[0][0]+distance*0.1+2*h};
  physics::max_boundary = {range_boundaries[1][0]-distance*0.1-2*h};

  clog_one(info) << "Boundaries=" << physics::min_boundary << ";"
         << physics::max_boundary << std::endl;

  bs.update_iteration();

#ifdef OUTPUT
  bs.write_bodies("output_sodtube",physics::iteration);
#endif

  ++physics::iteration; 
  do
  { 
    analysis::screen_output();
    MPI_Barrier(MPI_COMM_WORLD);


    // Compute and prepare the tree for this iteration 
    // - Compute the Max smoothing length 
    // - Compute the range of the system using the smoothinglength
    // - Cmopute the keys 
    // - Distributed qsort and sharing 
    // - Generate and feed the tree
    // - Exchange branches for smoothing length 
    // - Compute and exchange ghosts in real smoothing length 
    bs.update_iteration();

    clog_one(trace) << "compute_density_pressure_soundspeed" << std::flush; 
    bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
    clog_one(trace) << ".done" << std::endl;

    bs.update_neighbors();
   
    clog_one(trace) << "Hydro acceleration" << std::flush; 
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    clog_one(trace) << ".done" << std::endl;
 
    clog_one(trace) << "Internalenergy" << std::flush; 
    bs.apply_in_smoothinglength(physics::compute_dudt);
    clog_one(trace) << ".done" << std::endl;
   
    if(physics::iteration==1){ 
      clog_one(trace) << "leapfrog" << std::flush; 
      bs.apply_all(physics::leapfrog_integration_first_step);
      clog_one(trace) << ".done" << std::endl;
    }else{
      clog_one(trace) << "leapfrog" << std::flush; 
      bs.apply_all(physics::leapfrog_integration);
      clog_one(trace) << ".done" << std::endl;
    }

    clog_one(trace) << "dudt integration" << std::flush; 
    bs.apply_all(physics::dudt_integration);
    clog_one(trace) << ".done" << std::endl;

#ifdef OUTPUT_ANALYSIS
    // Compute the analysis values based on physics 
    bs.get_all(analysis::compute_lin_momentum);
    bs.get_all(analysis::compute_total_mass);
    // Only add the header in the first iteration
    analysis::scalar_output("scalar_sodtube.dat");
#endif

#ifdef OUTPUT
    if(physics::iteration % iteroutput == 0){ 
      bs.write_bodies("output_sodtube",physics::iteration/iteroutput);
    }
#endif
    ++physics::iteration;
    physics::totaltime += physics::dt;
    
  }while(physics::iteration<numberiterations);
}

flecsi_register_mpi_task(mpi_init_task);

void usage()
{
  clog_one(warn)<<"./sodtube [number of iterations]"<<std::endl;
}

void 
specialization_tlt_init(int argc, char * argv[]){
  
  // Default start at iteration 0
  usage();

  int numberiterations = 100;
  if(argc == 2){
    numberiterations = atoi(argv[1]);
  }
  clog_one(trace) << "In user specialization_driver" << std::endl;
  flecsi_execute_mpi_task(mpi_init_task,numberiterations); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  clog_one(trace) << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


