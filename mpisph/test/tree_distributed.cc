#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>

#include "io.h"
#include "tree_colorer.h"
//#include "tree_geometry.h"

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

using geometry_t = flecsi::topology::tree_geometry<double, 3>;

TEST(tree_distribution, distribution) {

  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);

  // Perform the computation on one process to compute the neighbors
  std::vector<body> bodies_check;
  int64_t totalnbodies_check,  localnbodies_check;
  io::inputDataHDF5(bodies_check,"tree_distributed_data.h5part","",
      totalnbodies_check,localnbodies_check,0,MPI_COMM_SELF);
  range_t range_check;
  range_check[0] = bodies_check[0].coordinates();
  range_check[1] = bodies_check[0].coordinates();
  // Compute the range
  for(size_t i = 0; i < bodies_check.size(); ++i){
    for(size_t d = 0; d < gdimension; ++d){

      if(bodies_check[i].coordinates()[d]+
          bodies_check[i].radius()>range_check[1][d])
        range_check[1][d] = bodies_check[i].coordinates()[d]+
                      bodies_check[i].radius();

      if(bodies_check[i].coordinates()[d]-
          bodies_check[i].radius()<range_check[0][d])
        range_check[0][d] = bodies_check[i].coordinates()[d]-
                      bodies_check[i].radius();
    }
  }

  // Compute the keys
  for(size_t i = 0; i < bodies_check.size(); ++i){
    bodies_check[i].set_key(entity_key_t(range_check,bodies_check[i].coordinates()));
  }

  // Sort
  #ifdef BOOST_PARALLEL
      boost::sort::block_indirect_sort(
  #else
      std::sort(
  #endif
          bodies_check.begin(),bodies_check.end(),
          [](auto& left, auto& right){
          if(left.key() < right.key()){
            return true;
          }
          if(left.key() == right.key()){
            return left.id() < right.id();
          }
          return false;
        }); // sort

  // Compute the neighbors in N^2 way
  for(auto& b1: bodies_check)
  {
    int neighbors = 0;
    for(auto& b2: bodies_check)
    {
      if(geometry_t::within_square(
            b1.coordinates(),
            b2.coordinates(),
            b1.radius(),b2.radius()))
      {
        ++neighbors;
      }
    }
    b1.set_neighbors(neighbors);
  }

  assert(localnbodies_check == totalnbodies_check);

  tree_colorer<double,gdimension> tc;
  tree_topology_t t;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  rank || clog(trace)<<"#proc: "<<size<<std::endl;

  int64_t totalnbodies;
  int64_t localnbodies;
  rank || clog(trace)<<"Reading particles"<<std::endl;
  io::inputDataHDF5(t.entities(),"tree_distributed_data.h5part","",
      totalnbodies,localnbodies,0);

  // Check the total number of bodies
  int64_t checknparticles = t.entities().size();
  MPI_Allreduce(MPI_IN_PLACE,&checknparticles,1,MPI_INT64_T,
    MPI_SUM,MPI_COMM_WORLD);
  assert(checknparticles==totalnbodies);

  // Compute the range
  rank || clog(trace)<<"Computing range"<<std::endl;
  range_t range;
  tc.mpi_compute_range(t.entities(),range);
  rank || clog(trace)<<"range: "<<range[0]<<";"<<range[1]<<std::endl;

  // Set the tree range
  rank || clog(trace)<<"Setting tree range"<<std::endl;
  t.set_range(range);
  // Compute the keys
  rank || clog(trace)<<"Compute keys"<<std::endl;
  t.compute_keys();
  // 2. Distribution
  rank || clog(trace)<<"Distributing "<<totalnbodies<<" particles"<<std::endl;
  tc.mpi_qsort(t.entities(),totalnbodies);

  rank || clog(trace)<<"Building tree "<<std::endl;

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

  rank || clog(trace)<<"Sharing edge"<<std::endl;
  // 4. Share the edge particles and compute COFM
  t.share_edge();
  t.cofm(t.root(),0.,false);

  rank || clog(trace)<<"Sharing branches"<<std::endl;
  // 5. Exchange the branches and recompute cofm
  std::vector<range_t> rangeposproc;
  tc.mpi_branches_exchange(t,t.entities(),rangeposproc,range);
  t.cofm(t.root(),0.,false);

  // Check the range of all the processes
  std::vector<range_t> root_ranges(size);
  range_t root_range = {t.root()->bmin(),t.root()->bmax()};
  MPI_Allgather(&root_range,sizeof(range_t),MPI_BYTE,
    &(root_ranges[0]),sizeof(range_t),MPI_BYTE,MPI_COMM_WORLD);

  int rk = 0 ;
  for(auto r: root_ranges)
  {
    assert(r[0] == root_range[0]);
    assert(r[1] == root_range[1]);
  }

  // Reduction on the local particles
  std::vector<int64_t> nparticles(size);
  MPI_Allgather(&localnbodies,1,MPI_INT64_T,&(nparticles[0]),1,MPI_INT64_T,
    MPI_COMM_WORLD);
  // Prefix scan
  std::vector<int64_t> nparticles_offset(size);
  std::partial_sum(nparticles.begin(),nparticles.end(),&nparticles_offset[0]);
  nparticles_offset.insert(nparticles_offset.begin(),0);

  // 6. Perform a several traversal and compute the neighbors
  int64_t ncritical = 1;
  bool variable_sph = true;
  t.apply_sub_cells(t.root(),0.,ncritical,variable_sph,
    [](body_holder * srch,std::vector<body_holder*>& nbh)
    {
      srch->getBody()->set_neighbors(nbh.size());
    });

  // Compare the neighbors with the result on one process
  // Find my start in the global particles for comparison
  size_t i = nparticles_offset[rank];
  for(auto bi: t.entities())
  {
    assert(bi.id() == bodies_check[i].id());
    if(bi.neighbors() != bodies_check[i].neighbors())
    {
      std::cerr<<bi<<" N = "<<bi.neighbors()<<std::endl;
      std::cerr<<bodies_check[i]<<" N = "<<bodies_check[i].neighbors()<<std::endl;
    }
    assert(bi.neighbors() == bodies_check[i].neighbors());
    ++i;
  }

  MPI_Finalize();
}
