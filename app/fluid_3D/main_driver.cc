/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*\

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
#include "fluids_physics.h"
#include "analysis.h"

namespace flecsi{
namespace execution{

void
mpi_init_task(int startiteration){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int noutput = startiteration+1;

  body_system<double,gdimension> bs;
  
  //bs.read_bodies("hdf5_fluid_3D.h5part",bs,startiteration);
  
  // Use the user reader defined in the physics file 
  bs.read_bodies("hdf5_fluid_3D.h5part",startiteration);
  physics::init_physics("hdf5_fluid_3D.h5part");

  physics::maxtime = 4.0;
  remove("output_fluid_3D.h5part");

  bs.update_iteration();

  // For OMP time, just used on 0, initialized for others
  double start = 0;

  // Apply EOS once 
  MPI_Barrier(MPI_COMM_WORLD);
  rank|| clog(trace)<<"Init EOS"<<std::flush; 
  start = omp_get_wtime(); 
  bs.apply_all(physics::EOS); 
  MPI_Barrier(MPI_COMM_WORLD);
  rank|| clog(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

  // Apply EOS once 
  rank|| clog(trace)<<"Init Density Velocity"<<std::flush; 
  start = omp_get_wtime(); 
  bs.apply_all(physics::init_density_velocity); 
  MPI_Barrier(MPI_COMM_WORLD);
  rank|| clog(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

#ifdef OUTPUT
  bs.write_bodies("output_fluid_3D",startiteration);
#endif

  ++physics::iteration; 
  do
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    analysis::screen_output(rank);

    // Compute and prepare the tree for this iteration 
    // - Compute the Max smoothing length 
    // - Compute the range of the system using the smoothinglength
    // - Cmopute the keys 
    // - Distributed qsort and sharing 
    // - Generate and feed the tree
    // - Exchange branches for smoothing length 
    // - Compute and exchange ghosts in real smoothing length 
    bs.update_iteration();
    
    // Do the Sod Tube physics
    rank|| clog(trace)<<"Accel"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_in_smoothinglength(physics::compute_accel);
    rank|| clog(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;


    rank|| clog(trace)<<"DT"<<std::flush; 
    start = omp_get_wtime(); 
    bs.get_all(physics::update_dt); 
    rank|| clog(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    // Reduce the local DT 
    MPI_Allreduce(MPI_IN_PLACE,&physics::dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
    rank|| clog(trace)<<"New DT = "<<physics::dt<<std::endl;

    rank|| clog(trace)<<"Verlet integration EOS"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_all(physics::verlet_integration_EOS);
    rank|| clog(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    physics::totaltime += physics::dt;
    rank|| clog(trace)<<"Total time="<<physics::totaltime<<"s / "<<physics::maxtime<<"s"<<std::endl;

#ifdef OUTPUT
    if(physics::totaltime > noutput*physics::outputtime){
      bs.write_bodies("output_fluid_3D",noutput);
      noutput++;
    }
#endif
    ++physics::iteration;
    
  }while(physics::iteration<40);
}

flecsi_register_mpi_task(mpi_init_task, flecsi::execution);

void 
specialization_tlt_init(int argc, char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  // Default start at iteration 0
  int startiteration = 0;
  if(argc == 2){
    startiteration = atoi(argv[1]);
  }

  rank|| clog(trace)<< "In user specialization_driver" << std::endl;

  flecsi_execute_mpi_task(mpi_init_task,flecsi::execution, startiteration); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  rank|| clog(trace) << "In user driver" << std::endl;
} // driver


} // namespace execution
} // namespace flecsi


