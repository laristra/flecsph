/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>
//#include <H5hut.h>

#include "user.h"
#include "sodtube.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"
#include "io.h"

using namespace io;
//
// help message
//
void print_usage() {
  clog(warn)
      << "Change the velocity for relaxed simulation"
      << gdimension << "D" << std::endl
      << "Usage: ./RT_XD_velocity <input_file_prefix> <output_file_prefix> "
      << "<iteration>" << std::endl;
}


//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace param;

  // check options list: exactly three options are allowed
  if (argc != 4) {
    print_usage();
    exit(0);
  }

  // launch MPI
  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // anything other than 2D is not implemented yet
  assert (gdimension == 2);
  assert (domain_type == 0);

  // set simulation parameters
  //param::mpi_read_params(argv[1]);
  //set_derived_params();

  // Input the bodies
  // Input data fro HDF5 File
  std::vector<std::pair<entity_key_t,body>> bodies;
  int64_t totalnbodies, nbodies;

  io::inputDataHDF5(bodies,argv[1],argv[2],totalnbodies,nbodies,
   atoi(argv[3]));

  for(auto& b: bodies)
  {
    point_t vel = b.second.getVelocity();
    point_t pos = b.second.coordinates();
    // Add velocity perturbation a-la Price (2008)
    vel[1] = 0.01*(1 + cos(4*M_PI*pos[0]))*(1 + cos(3*M_PI*pos[1]))/4.;
    b.second.setVelocity(vel);
  } // for part=0..nparticles

  io::outputDataHDF5(bodies,argv[2],0,0.);

  MPI_Finalize();
  return 0;
}
