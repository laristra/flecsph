#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>

#include "tree_colorer.h"

// Number of particles total
#define N 400
#define RMINX 0.
#define RMINY 0.
#define RMINZ 0.
#define RMAXX 0.1
#define RMAXY 0.1
#define RMAXZ 0.1
#define HMAX 0.001
#define HMIN 0.0001

using namespace std;
using namespace flecsi;
using namespace topology;

std::ostream&
operator<<(
  std::ostream& ostr,
  const entity_key_t id)
{
  id.output_(ostr);
  return ostr;
}

double uniform(){
  return double(rand())/RAND_MAX;
}

double uniform(double a, double b){
  return a + (b - a) * uniform();
}

namespace flecsi{
namespace execution{
  void driver(int argc, char* argv[]){
  }
}
}

TEST(tree_distribution, distribution) {
  tree_colorer<double,gdimension> tc;
  tree_topology_t t;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  srand(rank);

  rank || clog(trace)<<"#proc: "<<size<<std::endl;

  size_t totalparticles = N;
  size_t localparticles = N / size;
  if(rank == size - 1 ){
    localparticles = totalparticles - (N/size)*(size-1);
  }
  rank || clog(trace)<<"#part: "<<localparticles<<" ("<<
    totalparticles - (N/size)*(size-1)<<")"<<std::endl;

  t.entities().resize(localparticles);

  double mass = 1.0;
  range_t range = {point_t{RMINX,RMINY,RMINZ},point_t{RMAXX,RMAXY,RMAXZ}};
  rank || clog(trace)<<"Range: "<<range[0]<<"-"<<range[1]<<std::endl;

  int64_t first_id = N/size*rank+1;
  // 1. generate particles
  for(size_t i = 0; i < localparticles; ++i){
    point_t p = {uniform(RMINX, RMAXX), uniform(RMINY, RMAXY),
        uniform(RMINZ, RMAXZ)};
    double radius = uniform(HMIN,HMAX);
    t.entities()[i].set_coordinates(p);
    t.entities()[i].set_mass(mass);
    t.entities()[i].set_id(i+first_id);
    t.entities()[i].set_radius(radius);
    t.entities()[i].setDensity(0);
    t.entities()[i].setPressure(0);
  }
  // Set the tree range
  rank || clog(trace)<<"Setting tree range"<<std::endl;
  t.set_range(range);
  // Compute the keys
  rank || clog(trace)<<"Compute keys"<<std::endl;
  t.compute_keys();
  // 2. Distribution
  rank || clog(trace)<<"Distributing particles"<<std::endl;
  tc.mpi_qsort(t.entities(),totalparticles);

  // Check the distribution
  rank || clog(trace) << "Done"<<std::endl<<std::flush;
}


TEST(tree_distribution, create_tree) {

  tree_colorer<double,gdimension> tc;
  tree_topology_t t;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  srand(rank);

  rank || clog(trace)<<"#proc: "<<size<<std::endl;

  size_t totalparticles = N;
  size_t localparticles = N / size;
  if(rank == size - 1 ){
    localparticles = totalparticles - (N/size)*(size-1);
  }

  t.entities().resize(localparticles);

  double mass = 1.0;
  range_t range = {point_t{RMINX,RMINY,RMINZ},point_t{RMAXX,RMAXY,RMAXZ}};
  rank || clog(trace)<<"Range: "<<range[0]<<"-"<<range[1]<<std::endl;

  int64_t first_id = N/size*rank+1;
  // 1. generate particles
  for(size_t i = 0; i < localparticles; ++i){
    point_t p = {uniform(RMINX, RMAXX), uniform(RMINY, RMAXY),
        uniform(RMINZ, RMAXZ)};
    double radius = uniform(HMIN,HMAX);
    t.entities()[i].set_coordinates(p);
    t.entities()[i].set_mass(mass);
    t.entities()[i].set_id(i+first_id);
    t.entities()[i].set_radius(radius);
    t.entities()[i].setDensity(0);
    t.entities()[i].setPressure(0);
  }
  // Set the tree range
  t.set_range(range);
  // Compute the keys
  t.compute_keys();
  // 2. Distribution
  tc.mpi_qsort(t.entities(),totalparticles);

  // 3. Compute the local tree
  for(auto& bi:  t.entities()){
    bi.set_owner(rank);
    auto id = t.make_entity(bi.key(),bi.coordinates(),
      &(bi),rank,bi.mass(),bi.id(),bi.radius());
    t.insert(id);
    auto nbi = t.get(id);

    assert(nbi->global_id() == bi.id());
    assert(nbi->getBody() != nullptr);
    assert(nbi->is_local());
  }
}
