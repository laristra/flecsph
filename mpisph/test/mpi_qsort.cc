
#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "tree_colorer.h"

using namespace ::testing;


// Rule for body equals == same position
inline 
bool 
operator== (
  const body& b1, 
  const body& b2)
{
  return b1.getPosition()==b2.getPosition();    
};

namespace flecsi{
  namespace execution{
    void driver(int argc, char* argv[]){
    }
  }
}

TEST(tree_colorer, mpi_qsort){
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  srand(time(NULL)*rank);
  tree_colorer<double,gdimension> tc;

  // Generating the particles randomly on each process
  size_t nparticles = 10000;
  size_t nparticlesperproc = nparticles/size; 
  double maxbound = 1.0; // Particles positions between [0,1]
  // Adjust for last one 
  if(rank == size-1){
    nparticlesperproc = (nparticles-nparticlesperproc*(size-1));
  }
  clog_one(info)<<"Generating "<<nparticles<<std::endl;

  std::cout<<"Rank "<<rank<<": "<<nparticlesperproc<<" particles"<<std::endl;

  // Range to compute the keys 
  std::array<point_t,2> range;
  range[0] = point_t{};
  range[1] = point_t{maxbound,maxbound,maxbound};
  entity_key_t::set_range(range); 
  std::vector<std::pair<entity_key_t,body>> bodies(nparticlesperproc);
  // Create the bodies and keys 
  for(size_t i=0;i<nparticlesperproc;++i){
    // Random x, y and z
    bodies[i].second.setPosition(
      point_t{(double)rand()/(double)RAND_MAX*(maxbound),
      (double)rand()/(double)RAND_MAX*(maxbound),
      (double)rand()/(double)RAND_MAX*(maxbound)}
    );
 
    // Compute the key 
    bodies[i].first = entity_key_t(/*range,*/bodies[i].second.getPosition());
  }

  // Gather all the particles everywhere and sort locally 
  std::vector<std::pair<entity_key_t,body>> checking(nparticles);
  MPI_Allgather(
    &bodies[0],
    nparticlesperproc*sizeof(std::pair<entity_key_t,body>),
    MPI_BYTE,
    &checking[0],
    nparticlesperproc*sizeof(std::pair<entity_key_t,body>),
    MPI_BYTE,
    MPI_COMM_WORLD);

  // Sort it locally base on the keys 
  std::sort(checking.begin(),checking.end(),
      [](auto& left, auto& right){
        return left.first < right.first;
      });

  // Extract the subset of this process 
  std::vector<std::pair<entity_key_t,body>> my_checking(
    checking.begin()+rank*(nparticles/size),
    checking.begin()+rank*nparticlesperproc+nparticlesperproc);

  // Use the mpi_qsort
  tc.mpi_qsort(bodies,nparticles); 
	
	// Compare the results with all processes particles subset 
  ASSERT_TRUE(my_checking == bodies);
}
