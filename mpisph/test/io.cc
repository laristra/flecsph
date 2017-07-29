
#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "tree_colorer.h"
#include "io.h"

using namespace ::testing;

namespace flecsi{
  namespace execution{
    void driver(int argc, char* argv[]){
    }
  }
}

TEST(io, open_file){
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  MPI_Barrier(MPI_COMM_WORLD);
  std::cout<<rank<<"/"<<size<<std::endl;

  // Open the file
  auto dataFile = H5OpenFile(
      "io_test.h5part",
      H5_O_RDONLY|H5_VFD_INDEPENDENT, // Flag to not use mpiposix
      MPI_COMM_WORLD);

  // Close the file 
  H5CloseFile(dataFile);
}

TEST(io, read_file){
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


}

TEST(io, write_file){ 
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

}
