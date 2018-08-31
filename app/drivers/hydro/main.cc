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
 * @file main.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Main function, start MPI with Gasnet. Then launch fleCSI runtime.  
 */

#include <flecsi.h>
#include "flecsi/execution/execution.h"
#include "flecsi/concurrency/thread_pool.h"

#if FLECSI_RUNTIME_MODEL == FLECSI_RUNTIME_MODEL_legion
  #include <mpi.h>
  #include <legion.h>
#endif

int main(int argc, char * argv[]){
  
  int provided;
 
  // Normal way 
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (provided < MPI_THREAD_MULTIPLE)
    printf("ERROR: Your implementation of MPI does not support "
     "MPI_THREAD_MULTIPLE which is required for use of the "
     "GASNet MPI conduit with the Legion-MPI Interop!\n");
  assert(provided == MPI_THREAD_MULTIPLE);

  auto retval = flecsi::execution::context_t::instance().initialize(argc,argv);

  return retval;

}


