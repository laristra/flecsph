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
  
  clog(trace)<<"Creating bodies"<<std::endl<<std::flush;
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

  clog(trace)<<"Computing keys"<<std::endl<<std::flush;
  // Compute the keys 
  for(auto& bi:  bodies){
    bi.first = entity_key_t(tree->range(),bi.second.coordinates());
  }
 
  // Bodies holders
  std::vector<body_holder*> bodies_holders; 
 
 // Set the size of the vector 
 tree->set_entities_vector_size(nbodies); 
  
 clog(trace)<<"Adding in tree"<<std::endl<<std::flush; 
  // Add my local bodies in my tree 
  for(auto& bi:  bodies){
    auto nbi = tree->make_entity(bi.second.getPosition(),&(bi.second),0,
      bi.second.getMass(),bi.second.getId(),bi.second.getSmoothinglength());
    tree->insert(nbi); 
    bodies_holders.push_back(nbi);
    assert(nbi->global_id() == bi.second.id());
    assert(nbi->getBody() != nullptr);
    assert(nbi->is_local());
  }

  // Destroy the tree 
  delete tree; 
}
