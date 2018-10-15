#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>

#include "tree.h"

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

TEST(tree_topology, neighbors_sphere) {
  tree_topology_t t;

  vector<body_holder*> ents;

  size_t n = 5000;
  double mass = 1.0;

  for(size_t i = 0; i < n; ++i){
    point_t p = {uniform(0, 1), uniform(0, 1), uniform(0, 1)};
    auto e = t.make_entity(p,nullptr,0,mass,0);
    t.insert(e);
    ents.push_back(e);
  }

  t.update_branches(0.1);

  ASSERT_TRUE(t.root()->mass() == n*mass);

  for(size_t i = 0; i < n; ++i){
    auto ent = ents[i];

    auto ns = t.find_in_radius_b(ent->coordinates(), 0.10);

    set<body_holder*> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder*> s2;

    for(size_t j = 0; j < n; ++j){
      auto ej = ents[j];

      if(distance(ent->coordinates(), ej->coordinates()) < 0.10){
        s2.insert(ej);
      }
    }

    //std::cout<<s1.size()<<" && " << s2.size()<<std::endl;

    ASSERT_TRUE(s1 == s2);
  }
}

TEST(tree_topology, neighbors_box) {
  tree_topology_t t;

  vector<body_holder*> ents;

  size_t n = 5000;
  double mass = 1.0;
  
  point_t max;
  point_t min;

  for(size_t i = 0; i < n; ++i){
    point_t p = {uniform(0, 1), uniform(0, 1), uniform(0, 1)};
    auto e = t.make_entity(p,nullptr,0,mass,0);
    t.insert(e);
    ents.push_back(e);
  }

  t.update_branches(0.1);
  
  ASSERT_TRUE(t.root()->mass() == n*mass);

  for(size_t i = 0; i < n; ++i){
    auto ent = ents[i];
	
	for(size_t d = 0; d < gdimension; ++d ){
		max[d] = ent->coordinates()[d]+0.1;
		min[d] = ent->coordinates()[d]-0.1;
	}
    auto ns = t.find_in_box_b(min,max);

	//std::cout<<ns.size()<<std::endl;
	
    set<body_holder*> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder*> s2;

    for(size_t j = 0; j < n; ++j){
      auto ej = ents[j];

	  bool in_box = true;
	  for(size_t d = 0; d < gdimension; ++d ){
		if(ej->coordinates()[d] > max[d] || ej->coordinates()[d] < min[d]){
			in_box = false; 
			break;
		}		
	  }
	  
      if(in_box){
        s2.insert(ej);
      }
    }

    //std::cout<<s1.size()<<" && " << s2.size()<<std::endl;

    ASSERT_TRUE(s1 == s2);
  }
}

TEST(tree_topology,smoothing){
  size_t niter = 20; 
  size_t nparticles_line = 4;
  size_t nparticles = nparticles_line*nparticles_line*nparticles_line;
  double distance = 0.1;
  double epsilon = 0.000001;

  double h = distance/2.-10*epsilon;

  //std::cout<<"Distance="<<distance<<" H="<<h<<" Epsilon="<<epsilon<<std::endl;

  // Generate the particles 
  vector<body_holder*> ents;
  double line =  0.;
  double col  =  0.;
  double depth = 0.;
  double mass = 1.0;
  point_t min = {0.-distance,0.-distance,0.-distance};
  point_t max = { (nparticles_line-1)*distance+distance,
                  (nparticles_line-1)*distance+distance,
                  (nparticles_line-1)*distance+distance};
  //std::cout<<"range="<<min<<max<<std::endl;
  tree_topology_t t(min,max);
  for(size_t part_line=0; part_line<nparticles_line;++part_line){
      col = 0.;
    for(size_t part_col=0; part_col<nparticles_line;++part_col){
      depth = 0.;
      for(size_t part_depth=0; part_depth<nparticles_line;++part_depth){
        point_t position = {line,col,depth};
        auto e = t.make_entity(position,nullptr,0,mass,0);
        t.insert(e);
        ents.push_back(e);
        depth += distance; 
        //std::cout<<*e<<std::endl;
      }
      col += distance; 
    }
    line += distance; 
  }

  t.update_branches(2*h);
  ASSERT_TRUE(t.root()->mass() == nparticles*mass); 


  for(size_t iter=0;iter < niter ; ++iter){
    //std::cout<< iter <<" h="<<h<<std::endl;
    // For each particle search in radius 
    for(size_t i = 0; i < nparticles; ++i){
      auto ent = ents[i];
      auto ns = t.find_in_radius_b(ent->coordinates(), 2*h);
      auto vec = ns.to_vec();

      vector<body_holder*> s1;
      s1.insert(s1.begin(),vec.begin(), vec.end());

      vector<body_holder*> s2;

      for(size_t j = 0; j < nparticles; ++j){
        auto ej = ents[j];

        if(flecsi::distance(ent->coordinates(), ej->coordinates()) <= 2*h){
          s2.push_back(ej);
        }
      }

      ASSERT_TRUE(s1.size() > 0.);
      ASSERT_TRUE(s2.size() > 0.);

      // Sort the vectors 
      std::sort(s1.begin(),s1.end());
      std::sort(s2.begin(),s2.end());  
      ASSERT_TRUE(s1 == s2);
    }
    // Change h
    h += epsilon; 
  }
}

#if 0
TEST(tree_topology,same_key){
  size_t nparticles_line = 20;
  size_t nparticles = nparticles_line*nparticles_line*nparticles_line;
  double distance = 1.0e-15;
  double space = 1.0;
  double h = distance*2.;

  // Generate the particles 
  vector<body_holder*> ents;
  vector<entity_key_t> keys; 
  double line =  0.;
  double col  =  0.;
  double depth = 0.;
  double mass = 1.0;
  point_t min = {-space,-space,-space};
  point_t max = {space,space,space};
  std::array<point_t,2> range = {min,max}; 
  tree_topology_t t(min,max);
  entity_key_t::set_range(range); 
  for(size_t part_line=0; part_line<nparticles_line;++part_line){
      col = 0.;
    for(size_t part_col=0; part_col<nparticles_line;++part_col){
      depth = 0.;
      for(size_t part_depth=0; part_depth<nparticles_line;++part_depth){
        point_t position = {line,col,depth};
        auto e = t.make_entity(position,nullptr,0,mass,0);
        t.insert(e);
        ents.push_back(e);
        entity_key_t tmp = entity_key_t(/*range,*/position);
        keys.push_back(tmp);
        depth += distance; 
        //std::cout<<position<<" key="<<tmp<<std::endl;
      }
      col += distance; 
    }
    line += distance; 
  }

  std::sort(keys.begin(),keys.end());
  if(!(keys.end() == std::unique(keys.begin(),keys.end()))){
    std::cout<<"Key colision"<<std::endl;
  } 

  t.update_branches(2*h); 
  ASSERT_TRUE(t.root()->getMass() == nparticles*mass); 


  for(size_t i = 0; i < nparticles; ++i){
    auto ent = ents[i];
    auto ns = t.find_in_radius_b(ent->coordinates(), 2*h);
    auto vec = ns.to_vec();

    vector<body_holder*> s1;
    s1.insert(s1.begin(),vec.begin(), vec.end());

    vector<body_holder*> s2;

    for(size_t j = 0; j < nparticles; ++j){
      auto ej = ents[j];

      if(flecsi::distance(ent->coordinates(), ej->coordinates()) <= 2*h){
        s2.push_back(ej);
      }
    }

    ASSERT_TRUE(s1.size() > 0.);
    ASSERT_TRUE(s2.size() > 0.);

    // Sort the vectors 
    std::sort(s1.begin(),s1.end());
    std::sort(s2.begin(),s2.end());  
    ASSERT_TRUE(s1 == s2);
  } 
}
#endif 
