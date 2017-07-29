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

#include "physics/eos_analytics.h"
//#include "physics.h"

namespace flecsi{
namespace execution{

void
mpi_init_task(){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 
  std::cout<<rank<<"/"<<size<<" Opening"<<std::endl<<std::flush;
  
  Flecsi_Sim_IO::HDF5ParticleIO test(
      "io.h5part",
      Flecsi_Sim_IO::READING,
      MPI_COMM_WORLD); 

  std::cout<<rank<<"/"<<size<<" Opening done, other version"
    <<std::endl<<std::flush;

  auto dataFile = H5OpenFile(
      "io.h5part", 
      H5_O_RDWR,
      MPI_COMM_WORLD);

  std::cout<<rank<<"/"<<size<<" Opening done"<<std::endl<<std::flush;
  std::cout<<rank<<"/"<<size<<" Closing"<<std::endl<<std::flush;

  H5CloseFile(dataFile);

  std::cout<<rank<<"/"<<size<<" Closing done"<<std::endl<<std::flush;

}

flecsi_register_mpi_task(mpi_init_task);

void 
specialization_tlt_init(int argc, char * argv[]){
  
  std::cout << "In user specialization_driver" << std::endl;
  flecsi_execute_mpi_task(mpi_init_task); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


