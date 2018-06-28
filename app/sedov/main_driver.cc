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

namespace flecsi{
namespace execution{

void
mpi_init_task(int totaliterations){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int totaliters = totaliterations;
  int iteroutput = 1;
  double totaltime = 0.0;
  double maxtime = 10.0;
  int iter = 0; 

  // Init if default values are not ok
  physics::dt = 0.001;
  physics::alpha = 1; 
  physics::beta = 2; 
  physics::do_boundaries = true;
  physics::stop_boundaries = true;
  physics::gamma = 5./3.;

  body_system<double,gdimension> bs;
  bs.read_bodies("hdf5_sedov.h5part",iter);

  double h = bs.getSmoothinglength();
  physics::epsilon = 0.01*h*h;

  auto range_boundaries = bs.getRange(); 
  point_t distance = range_boundaries[1]-range_boundaries[0];
  for(int i = 0; i < gdimension; ++i){
    distance[i] = fabs(distance[i]);
  }
  physics::min_boundary = {(0.1+2*h)*distance+range_boundaries[0]};
  physics::max_boundary = {-(0.1-2*h)*distance+range_boundaries[1]};
  clog(info) << "Limits: " << physics::min_boundary << " ; "
         << physics::max_boundary << std::endl;

  remove("output_sedov.h5part"); 

#ifdef OUTPUT
  bs.write_bodies("output_sedov",iter);
#endif

  double stopt, startt; 

  ++iter; 
  do
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    clog(info) << "#### Iteration " << iter << std::endl;
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
   
    clog(info) << "compute_density_pressure_soundspeed" << std::flush; 
    bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
    clog(info) << ".done" << std::endl;
    
    // Refresh the neighbors within the smoothing length 
    bs.update_neighbors(); 

    clog(info) << "Hydro acceleration" << std::flush; 
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    clog(info) << ".done" << std::endl;
 
    clog(info) << "Internalenergy"<<std::flush; 
    bs.apply_in_smoothinglength(physics::compute_dudt);
    clog(info) << ".done" << std::endl;
   
    if(iter==1){ 
      clog(info) << "leapfrog" << std::flush; 
      bs.apply_all(physics::leapfrog_integration_first_step);
      clog(info) << ".done" << std::endl;
    }else{
      if(rank==0)
      clog(info) << "leapfrog" << std::flush; 
      bs.apply_all(physics::leapfrog_integration);
      clog(info) << ".done" << std::endl;
    }

    clog(info) << "dudt integration" << std::flush; 
    bs.apply_all(physics::dudt_integration);
    clog(info) << ".done" << std::endl;

#ifdef OUTPUT
    MPI_Barrier(MPI_COMM_WORLD);
    startt = omp_get_wtime();
    if(iter % iteroutput == 0){ 
      bs.write_bodies("output_sedov",iter/iteroutput);
    }
    stopt = omp_get_wtime();
    clog(info) << "Output time: " << omp_get_wtime()-startt << "s" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    ++iter;
    
  }while(iter<totaliters);
}

flecsi_register_mpi_task(mpi_init_task);

void 
specialization_tlt_init(int argc, char * argv[]){
  
  // Default start at iteration 0
  int totaliterations = 100;
  if(argc == 2){
    totaliterations = atoi(argv[1]);
  }

  clog(trace) << "In user specialization_driver" << std::endl;

  flecsi_execute_mpi_task(mpi_init_task,totaliterations); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  clog(trace) << "In user driver" << std::endl;
} // driver


} // namespace execution
} // namespace flesci


