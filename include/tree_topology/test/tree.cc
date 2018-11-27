#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "../../tree.h"

using namespace flecsi;
using namespace topology;

namespace flecsi{
  namespace execution{
    void driver(int argc, char * argv[]){}
  }
}

const size_t dimension = gdimension;
using range_t = std::array<point_t,2>;

TEST(tree, add_entities){
  range_t range{point_t(0.,0.,0.),point_t(1.,1.,1.)};

  tree_topology_t * tree;
  tree = new tree_topology_t(range[0],range[1]);

  //clog(trace)<<"Creating bodies"<<std::endl<<std::flush;
  // Create bodies
  size_t nbodies = 1000;
  std::vector<std::pair<entity_key_t,body>> bodies(nbodies);
  size_t id = 0;
  for(size_t i{0}; i < nbodies; ++i){
    bodies[i].second.setPosition(point_t(
          (double)rand()/(double)RAND_MAX,
          (double)rand()/(double)RAND_MAX,
          (double)rand()/(double)RAND_MAX
        ));
    bodies[i].second.setMass((double)rand()/(double)RAND_MAX);
    bodies[i].second.setId(id++);
    bodies[i].second.setSmoothinglength((double)rand()/(double)RAND_MAX);
  }

  //clog(trace)<<"Computing keys"<<std::endl<<std::flush;
  // Compute the keys
  for(auto& bi:  bodies){
    bi.first = entity_key_t(tree->range(),bi.second.coordinates());
  }

 //clog(trace)<<"Adding in tree"<<std::endl<<std::flush;
  // Add my local bodies in my tree
  for(auto& bi:  bodies){
    //std::cout<<"Creating body"<<std::endl<<std::flush;
    auto id = tree->make_entity(bi.first,bi.second.getPosition(),&(bi.second),0,
      bi.second.getMass(),bi.second.getId(),bi.second.getSmoothinglength());
    //std::cout<<"Create ID:"<<id<<std::endl<<std::flush;
    tree->insert(id);
    //std::cout<<"Inserted"<<std::endl<<std::flush;
    auto nbi = tree->get(id);
    //std::cout<<"Got body:"<<*nbi<<std::endl<<std::flush;
    assert(nbi->global_id() == bi.second.id());
    assert(nbi->getBody() != nullptr);
    assert(nbi->is_local());
  }

  // Destroy the tree
  delete tree;
}
