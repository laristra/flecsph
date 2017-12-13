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
  
  int iter = startiteration; 
  int noutput = startiteration+1;
  int maxiter = 300;

  body_system<double,gdimension> bs;
  double maxtime = 1.0;
  double outputtime = 0.001;
  
  // Use the user reader defined in the physics file 
  bs.read_bodies("hdf5_bns_3D.h5part",startiteration);

  remove("output_bns_3D.h5part");

  bs.update_iteration();
    // For OMP time, just used on 0, initialized for others
  double start = 0;
  double start_iteration = 0;

  // Compute density, pressure, cs for next iteration
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    std::cout<<"compute_density_pressure_soundspeed"<<std::flush; 
    start = omp_get_wtime(); 
  }
  bs.apply_in_smoothinglength(
    physics::compute_density_pressure_soundspeed); 
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

#ifdef OUTPUT
  bs.write_bodies("output_bns_3D",startiteration);
#endif

  ++iter; 
  do
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    start_iteration = omp_get_wtime();
    if(rank==0)
      std::cout<<std::endl<<"#### Iteration "<<iter<<std::endl;
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
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Accel FMM"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.gravitation_fmm();
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    physics::dt = 0.001;

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Leapfrog integration"<<std::flush; 
      start = omp_get_wtime(); 
    }
    if(iter == 1){
      bs.apply_all(physics::leapfrog_integration_first_step);
    }else{
      bs.apply_all(physics::leapfrog_integration);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;


    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"compute_density_pressure_soundspeed"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_in_smoothinglength(
      physics::compute_density_pressure_soundspeed); 
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    physics::totaltime += physics::dt;

    if(rank == 0){
      std::cout<<"Total time="<<physics::totaltime<<"s / "<<maxtime<<"s"<<std::endl;
    }
#ifdef OUTPUT
    if(physics::totaltime > noutput*outputtime){
      bs.write_bodies("output_bns_3D",noutput);
      noutput++;
    }

    //if(iter % iteroutput == 0){ 
    //  bs.write_bodies("output_fluid",iter/iteroutput);
    //}
#endif
    ++iter;
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
      std::cout<<"Iteration time = "<<omp_get_wtime()-start_iteration<< "s" <<std::endl;
    }
  }while(iter<maxiter);
  //}while(physics::totaltime<physics::maxtime);
}

flecsi_register_mpi_task(mpi_init_task);

void 
specialization_tlt_init(int argc, char * argv[]){
  
  // Default start at iteration 0
  int startiteration = 0;
  if(argc == 2){
    startiteration = atoi(argv[1]);
  }

  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/
  /*std::string  filename(argv[1]);
  std::cout<<filename<<std::endl;*/
  flecsi_execute_mpi_task(mpi_init_task,startiteration); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


