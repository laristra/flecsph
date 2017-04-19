#ifndef tree_driver_h
#define tree_driver_h

#include <iostream>
#include <numeric> // For accumulate
#include <iostream>

#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

#include "physics.h"
#include "io.h"

namespace flecsi{
namespace execution{

void
mpi_task(/*const char * filename*/){
  const char * filename = "../data/data_binary_8338.txt";
  std::vector<body*> rbodies; // Body read by the process

  int rank; 
  int size; 
  int nbodies = 0;
  int totalnbodies = 0;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
  printf("Hi from process %d/%d, file %s\n",rank,size,filename);   
  // Read data from file, each process read a part of it 
  // For HDF5, no problems because we know the number of particles 
  // For txt format, work on the number of lines yet ... 
  io::inputDataTxtRange(rbodies,nbodies,totalnbodies,rank,size,filename); 

}

flecsi_register_task(mpi_task,mpi,index);

void 
specialization_driver(int argc, char * argv[]){
  if (argc!=2) {
    std::cerr << "Error not enough arguments\n"
        "Usage: tree <datafile>\n";
    exit(-1); 
  }


  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/

  flecsi_execute_task(mpi_task,mpi,index/*,filename*/); 
} // specialization driver

void 
driver(int argc, char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


#endif // tree_driver_h
