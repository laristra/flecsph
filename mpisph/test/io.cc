#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "utils.h"
#include "tree_colorer.h"
#include "io.h"


using namespace std;
using namespace flecsi;
using namespace topology;

namespace flecsi{
namespace execution{
  void driver(int argc, char* argv[]){
  }
}
}


TEST(io, write_N_read) {

  srand(time(NULL)); 
  int rank,size; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
  MPI_Comm_size(MPI_COMM_WORLD,&size); 

  const char * fileprefix = "io_utest"; 
  const char * filename = "io_utest.h5part";

  // Generate particles, write to file, read and compare 
  size_t n = 1000;
  double mass = 1.0;
  double u = 1.0;

  std::vector<std::pair<entity_key_t,body>> bodies(n); 

  for(size_t i = 0; i < n; ++i){
    point_t p = {
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX};
    // Create body to check the class too 
    bodies[i].second.setPosition(p);
    bodies[i].second.setMass(mass);
    bodies[i].second.setSmoothinglength(u);  
  }

  // Write that to file
  io::outputDataHDF5(bodies,fileprefix,0,0.); 

  std::vector<std::pair<entity_key_t,body>> rbodies; 
  int64_t totalnbodies = 0; 
  int64_t localnbodies = 0; 
  // Read this file 
  io::inputDataHDF5(rbodies,filename,totalnbodies,localnbodies,0);

  ASSERT_TRUE(localnbodies == n); 
  ASSERT_TRUE(totalnbodies == n*size); 

  // Remove the created file 
  remove(filename); 

}
