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

#include <bodies_system.h>

#include "default_physics.h"
#include "analysis.h"

#define OUTPUT_ANALYSIS

namespace flecsi{
namespace execution{

void
mpi_init_task(int totaliterations){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  clog_set_output_rank(0);
  
  int iteroutput = 10;
  //double totaltime = 0.0;
  double maxtime = 10.0;

  // Init if default values are not ok
  physics::dt = 0.001;
  physics::alpha = 2;   // Set to fit value in default_physics.h 
  physics::beta = 2; 
  physics::do_boundaries = false;
  physics::stop_boundaries = false;
  physics::gamma = 1.4; // Set to fit value in default_physics.h
  physics::epsilon = 0.01;

  body_system<double,gdimension> bs;
  bs.read_bodies("hdf5_sedov.h5part",physics::iteration);


  auto range_boundaries = bs.getRange(); 
  point_t distance = range_boundaries[1]-range_boundaries[0];
  for(int i = 0; i < gdimension; ++i){
    distance[i] = fabs(distance[i]);
  }
  double h = bs.getSmoothinglength();
  physics::min_boundary = {(0.1+2*h)*distance+range_boundaries[0]};
  physics::max_boundary = {-(0.1-2*h)*distance+range_boundaries[1]};
  clog_one(info) << "Limits: " << physics::min_boundary << " ; "
         << physics::max_boundary << std::endl;

  remove("output_sedov.h5part"); 

#ifdef OUTPUT
  bs.write_bodies("output_sedov",physics::iteration);
#endif

  double stopt, startt; 

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
    
    // Refresh the neighbors within the smoothing length 
    bs.update_neighbors(); 

    clog_one(trace) << "Hydro acceleration" << std::flush; 
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    clog_one(trace) << ".done" << std::endl;
 
    clog_one(trace) << "Internalenergy"<<std::flush; 
    bs.apply_in_smoothinglength(physics::compute_dudt);
    clog_one(trace) << ".done" << std::endl;
   
    if(physics::iteration==1){ 
      clog_one(trace) << "leapfrog" << std::flush; 
      bs.apply_all(physics::leapfrog_integration_first_step);
      clog_one(trace) << ".done" << std::endl;
    }else{
      if(rank==0)
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
    analysis::scalar_output("scalar_sedov.dat");
#endif

#ifdef OUTPUT
    MPI_Barrier(MPI_COMM_WORLD);
    startt = omp_get_wtime();
    if(physics::iteration % iteroutput == 0){ 
      bs.write_bodies("output_sedov",physics::iteration/iteroutput);
    }
    stopt = omp_get_wtime();
    clog_one(trace) << "Output time: " << omp_get_wtime()-startt << "s" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    ++physics::iteration;
    physics::totaltime += physics::dt;
    
  }while(physics::iteration<totaliterations);
}

flecsi_register_mpi_task(mpi_init_task);

void 
specialization_tlt_init(int argc, char * argv[]){
  
  // Default start at iteration 0
  int totaliterations = 100;
  if(argc == 2){
    totaliterations = atoi(argv[1]);
  }

  clog_one(trace) << "In user specialization_driver" << std::endl;

  flecsi_execute_mpi_task(mpi_init_task,totaliterations); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  clog_one(trace) << "In user driver" << std::endl;
} // driver


} // namespace execution
} // namespace flesci


