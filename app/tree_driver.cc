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

std::ostream&
operator<<(
  std::ostream& ostr,
  const entity_key_t& id
)
{
  id.output_(ostr);
  return ostr;
}

namespace flecsi{
namespace execution{

void
mpi_task(/*const char * filename*/){
  const char * filename = "../data/data_binary_8338.txt";
  //std::vector<body*> rbodies; // Body read by the process

  int rank; 
  int size; 
  int nbodies = 0;
  int totalnbodies = 0;
  tree_topology_t t;
  std::vector<std::pair<entity_key_t,body*>> rbodies;

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
  printf("%d/%d, file %s\n",rank,size,filename);   
  // Read data from file, each process read a part of it 
  // For HDF5, no problems because we know the number of particles 
  // For txt format, work on the number of lines yet ... 
  io::inputDataTxtRange(rbodies,nbodies,totalnbodies,rank,size,filename); 
  //std::cout << "Read done" << std::endl;

  // Compute the range to compute the keys 
  double maxx=-9999,maxy=-9999,maxz=-9999;
  double minx=9999,miny=9999,minz=9999;
  for(auto bi: rbodies){
    if(bi.second->coordinates()[0]>maxx)
      maxx = bi.second->coordinates()[0];
    if(bi.second->coordinates()[1]>maxy)
      maxy = bi.second->coordinates()[1];
    if(bi.second->coordinates()[2]>maxz)
      maxz = bi.second->coordinates()[2];
    if(bi.second->coordinates()[0]<minx)
      minx = bi.second->coordinates()[0];
    if(bi.second->coordinates()[1]<miny)
      miny = bi.second->coordinates()[1];
    if(bi.second->coordinates()[2]<minz)
      minz = bi.second->coordinates()[2];
  }
  // Do the MPI Reduction 
  MPI_Allreduce(MPI_IN_PLACE,&maxx,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
  MPI_Allreduce(MPI_IN_PLACE,&maxy,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
  MPI_Allreduce(MPI_IN_PLACE,&maxz,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
  MPI_Allreduce(MPI_IN_PLACE,&minx,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
  MPI_Allreduce(MPI_IN_PLACE,&miny,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
  MPI_Allreduce(MPI_IN_PLACE,&minz,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 

  point_t minposition = {minx,miny,minz};
  point_t maxposition = {maxx,maxy,maxz};

  std::cout << minposition << maxposition << std::endl;
  std::array<point_t,2> range = {minposition,maxposition};
  int maxdepth = 10;

  // The bodies are loaded
  // Compute the key and sort them 
  for(auto& bi: rbodies){
    bi.first = entity_key_t(range,bi.second->coordinates());
    //std::cout << bi.first << std::endl;
  }
  // Sort based on keys 
  std::sort(rbodies.begin(),rbodies.end(), [](auto &left,auto &right){
      return left.first < right.first; 
  });
  // Get the first and last key, send to everyone
  std::array<entity_key_t,2> locboundaries; 
  locboundaries[0] = rbodies.front().first;
  locboundaries[1] = rbodies.back().first; 
  std::cout <<rank<<"/"<<size<<" f: "<<locboundaries[0]<<" l: "
    <<locboundaries[1]<<std::endl;

  // Share the boundary keys, here we know it is a 64 unsigned int 
  // Make a vector for the right number of pairs 
  //vector<entity_key_t> globalboundaries;

  //MPI_Allgather(&locboundaries[0],2,);  

  // Create one vector per process with its keys 

  // Send the informations vio MPI_Alltoallv

  // Generate the local tree 

  // Do the research of ghost and shared 
  
  // Index everything

  // Register data and create the final tree ?
  MPI_Barrier(MPI_COMM_WORLD);

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
