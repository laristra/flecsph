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

#include "eos_analytics.h"
#include "physics.h"

namespace flecsi{
namespace execution{

void
mpi_init_task(int startiteration){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int totaliters = 100;
  int iteroutput = 1;
  double totaltime = 0.0;
  double maxtime = 10.0;
  int iter = startiteration; 

  // Init if default values are not ok
  physics::dt = 1.0e-10;
  physics::alpha = 1; 
  physics::beta = 2; 
  physics::stop_boundaries = true;
  physics::min_boundary = {0.1};
  physics::max_boundary = {1.0};
  physics::gamma = 5./3.;

  body_system<double,gdimension> bs;
  bs.read_bodies("bwd_id.h5part",startiteration);

  double h = bs.getSmoothinglength();
  physics::epsilon = 0.01*h*h;

#ifdef OUTPUT
  bs.write_bodies("output_bwd",iter);
#endif

  ++iter; 
  do
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    clog(info)<<"#### Iteration "<<iter<<std::endl;
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
   
    // Do the DWD physics
    clog(info)<<"Density"<<std::flush; 
    bs.apply_in_smoothinglength(physics::compute_density);
    clog(info)<<".done"<<std::endl;

    clog(info)<<"Pressure"<<std::flush; 
    bs.apply_all(physics::compute_pressure_wd);
    clog(info)<<".done"<<std::endl;

    clog(info)<<"Linear Momentum"<<std::flush; 
    //bs.apply_all(physics::compute_lin_momentum);
    clog(info)<<".done"<<std::endl;

    clog(info)<<"Soundspeed"<<std::flush; 
    bs.apply_all(physics::compute_soundspeed);
    clog(info)<<".done"<<std::endl;
    
    // Refresh the neighbors within the smoothing length 
    bs.update_neighbors(); 

    clog(info)<<"Hydro acceleration"<<std::flush; 
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    clog(info)<<".done"<<std::endl;
 
    clog(info)<<"Internalenergy"<<std::flush; 
    bs.apply_in_smoothinglength(physics::compute_dudt);
    clog(info)<<".done"<<std::endl; 
   
    if(iter==1){ 
      clog(info)<<"leapfrog"<<std::flush; 
      bs.apply_all(physics::leapfrog_integration_first_step);
      clog(info)<<".done"<<std::endl;
    }else{
      clog(info)<<"leapfrog"<<std::flush; 
      bs.apply_all(physics::leapfrog_integration);
      clog(info)<<".done"<<std::endl;
    }

    clog(info)<<"dudt integration"<<std::flush; 
    bs.apply_all(physics::dudt_integration);
    clog(info)<<".done"<<std::endl;

   
#ifdef OUTPUT
    if(iter % iteroutput == 0){ 
      bs.write_bodies("output_bwd",iter/iteroutput);
    }
#endif
    ++iter;
    
  }while(iter<totaliters);
}

flecsi_register_mpi_task(mpi_init_task);

void 
specialization_tlt_init(int argc, char * argv[]){
  
  // Default start at iteration 0
  int startiteration = 0;
  if(argc == 2){
    startiteration = atoi(argv[1]);
  }

  clog(warn) << "In user specialization_driver" << std::endl;

  flecsi_execute_mpi_task(mpi_init_task,startiteration); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  clog(warn) << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


