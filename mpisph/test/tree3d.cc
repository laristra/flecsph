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

TEST(tree_topology, neighbors_sphere_NORMAL) {
  tree_topology_t t;

  size_t n = 5000;
  double mass = 1.0;
  range_t range = {point_t{0,0,0},point_t{1,1,1}};

  for(size_t i = 0; i < n; ++i){
    point_t p = {uniform(0, 1), uniform(0, 1), uniform(0, 1)};
    auto e = t.make_entity(entity_key_t(range,p),p,nullptr,0,mass,0,0.1);
    t.insert(e);
  }


  t.post_order_traversal(t.root(),traversal_t::update_COM,
      0.00001,false);

  ASSERT_TRUE(t.root()->mass() == n*mass);

  for(size_t i = 0; i < n; ++i){
    auto ent = t.get(i);

    auto ns = t.find_in_radius(ent->coordinates(), 0.10,
        tree_geometry_t::within);

    set<body_holder*> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder*> s2;

    for(size_t j = 0; j < n; ++j){
      auto ej = t.get(j);

      if(distance(ent->coordinates(), ej->coordinates()) < 0.10){
        s2.insert(ej);
      }
    }

    ASSERT_TRUE(s1 == s2);
  }
}



TEST(tree_topology, neighbors_sphere_VARIABLE) {
  tree_topology_t t;

  size_t n = 5000;
  double mass = 1.0;
  range_t range = {point_t{0,0,0},point_t{1,1,1}};

  for(size_t i = 0; i < n; ++i){
    point_t p = {uniform(0, 1), uniform(0, 1), uniform(0, 1)};
    auto e = t.make_entity(entity_key_t(range,p),p,nullptr,0,mass,0,uniform(.1,.2));
    t.insert(e);
  }


  t.post_order_traversal(t.root(),traversal_t::update_COM,
      0.00001,false);

  ASSERT_TRUE(t.root()->mass() == n*mass);

  for(size_t i = 0; i < n; ++i){
    auto ent = t.get(i);

    auto ns = t.find_in_radius(ent->coordinates(), ent->h(),
        tree_geometry_t::within_square);

    set<body_holder*> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder*> s2;

    for(size_t j = 0; j < n; ++j){
      auto ej = t.get(j);
      double dist = distance(ent->coordinates(),ej->coordinates());
      if(dist*dist < (ent->h()+ej->h())*(ent->h()+ej->h())){
        s2.insert(ej);
      }
    }

    ASSERT_TRUE(s1 == s2);
  }
}

TEST(tree_topology, neighbors_box_NORMAL) {
  tree_topology_t t;

  size_t n = 5000;
  double mass = 1.0;
  range_t range = {point_t{0,0,0},point_t{1,1,1}};

  point_t max;
  point_t min;

  for(size_t i = 0; i < n; ++i){
    point_t p = {uniform(0, 1), uniform(0, 1), uniform(0, 1)};
    auto e = t.make_entity(entity_key_t(range,p),p,nullptr,0,mass,0,0.1);
    t.insert(e);
  }

  t.post_order_traversal(t.root(),traversal_t::update_COM,
      0.00001,false);

  ASSERT_TRUE(t.root()->mass() == n*mass);

  for(size_t i = 0; i < n; ++i){
    auto ent = t.get(i);

	for(size_t d = 0; d < gdimension; ++d ){
		max[d] = ent->coordinates()[d]+0.1;
		min[d] = ent->coordinates()[d]-0.1;
	}
    auto ns = t.find_in_box(min,max,tree_geometry_t::within_box);

    set<body_holder*> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder*> s2;

    for(size_t j = 0; j < n; ++j){
      auto ej = t.get(j);

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

    ASSERT_TRUE(s1 == s2);
  }
}


TEST(tree_topology, neighbors_box_VARIABLE) {
  tree_topology_t t;

  size_t n = 5000;
  double mass = 1.0;

  point_t max;
  point_t min;
  range_t range = {point_t{0,0,0},point_t{1,1,1}};

  for(size_t i = 0; i < n; ++i){
    point_t p = {uniform(0, 1), uniform(0, 1), uniform(0, 1)};
    auto e = t.make_entity(entity_key_t(range,p),p,nullptr,0,mass,0,uniform(.1,.2));
    t.insert(e);
  }

  t.post_order_traversal(t.root(),traversal_t::update_COM,
      0.00001,false);

  ASSERT_TRUE(t.root()->mass() == n*mass);

  for(size_t i = 0; i < n; ++i){
    auto ent = t.get(i);

	  for(size_t d = 0; d < gdimension; ++d ){
		  max[d] = ent->coordinates()[d]+0.00001+ent->h();
		  min[d] = ent->coordinates()[d]-0.00001-ent->h();
	  }
    auto ns = t.find_in_box(min,max,tree_geometry_t::intersects_sphere_box);

    set<body_holder*> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder*> s2;

    for(size_t j = 0; j < n; ++j){
      auto ej = t.get(j);

      if(tree_geometry_t::intersects_sphere_box(min,max,
        ej->coordinates(),ej->h())){
        s2.insert(ej);
      }
	  }
    ASSERT_TRUE(s1 == s2);
  }
}
