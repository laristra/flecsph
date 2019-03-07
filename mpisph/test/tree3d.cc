#include <cinchdevel.h>
#include <cinchtest.h>

#include <cmath>
#include <iostream>

#include "tree.h"

// Number of particles
#define N 2000
#define RMINX 0.
#define RMINY 0.
#define RMINZ 0.
#define RMAXX 0.01
#define RMAXY 0.01
#define RMAXZ 0.01
#define HMAX 0.001
#define HMIN 0.0001

using namespace std;
using namespace flecsi;
using namespace topology;

std::ostream &operator<<(std::ostream &ostr, const key_type id) {
  id.output_(ostr);
  return ostr;
}

double uniform() { return double(rand()) / RAND_MAX; }

double uniform(double a, double b) { return a + (b - a) * uniform(); }

namespace flecsi {
namespace execution {
void driver(int argc, char *argv[]) {}
} // namespace execution
} // namespace flecsi

TEST(tree_topology, neighbors_sphere_NORMAL) {
  tree_topology_t t;

  size_t n = N;
  double mass = 1.0;
  range_t range = {point_t{RMINX, RMINY, RMINZ}, point_t{RMAXX, RMAXY, RMAXZ}};
  std::cout << "Range: " << range[0] << "-" << range[1] << std::endl;

  for (size_t i = 0; i < n; ++i) {
    point_t p = {uniform(RMINX, RMAXX), uniform(RMINY, RMAXY),
                 uniform(RMINZ, RMAXZ)};
    auto e = t.make_entity(key_type(range, p), p, nullptr, 0, mass, 0, HMAX);
    t.insert(e);
  }

  std::cout << "Computing cofm";

  t.cofm(t.root(), 0, false);

  std::cout << ".done" << std::endl;

  ASSERT_TRUE(t.root()->mass() == n * mass);

  for (size_t i = 0; i < n; ++i) {
    auto ent = t.get(i);

    // std::cout<<"Entity "<<i+1<<"/"<<n<<" = "<<ent->key();
    auto ns =
        t.find_in_radius(ent->coordinates(), HMAX, tree_geometry_t::within);
    // std:cout<<" -> "<<ns.size()<<" nbrs.done"<<std::endl;

    set<body_holder *> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder *> s2;

    for (size_t j = 0; j < n; ++j) {
      auto ej = t.get(j);

      if (distance(ent->coordinates(), ej->coordinates()) < HMAX) {
        s2.insert(ej);
      }
    }

    ASSERT_TRUE(s1 == s2);
  }
}

TEST(tree_topology, neighbors_sphere_VARIABLE) {
  tree_topology_t t;

  size_t n = N;
  double mass = 1.0;
  range_t range = {point_t{RMINX, RMINY, RMINZ}, point_t{RMAXX, RMAXY, RMAXZ}};
  std::cout << "Range: " << range[0] << "-" << range[1] << std::endl;

  for (size_t i = 0; i < n; ++i) {
    point_t p = {uniform(RMINX, RMAXX), uniform(RMINY, RMAXY),
                 uniform(RMINZ, RMAXZ)};
    auto e = t.make_entity(key_type(range, p), p, nullptr, 0, mass, 0,
                           uniform(HMIN, HMAX));
    t.insert(e);
  }

  t.cofm(t.root(), 0, false);

  ASSERT_TRUE(t.root()->mass() == n * mass);

  for (size_t i = 0; i < n; ++i) {
    auto ent = t.get(i);

    auto ns = t.find_in_radius(ent->coordinates(), ent->radius(),
                               tree_geometry_t::within_square);

    set<body_holder *> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder *> s2;

    for (size_t j = 0; j < n; ++j) {
      auto ej = t.get(j);
      double dist = distance(ent->coordinates(), ej->coordinates());
      if (dist * dist < (ent->radius() + ej->radius()) *
                            (ent->radius() + ej->radius()) / 4.) {
        s2.insert(ej);
      }
    }

    ASSERT_TRUE(s1 == s2);
  }
}

TEST(tree_topology, neighbors_box_NORMAL) {
  tree_topology_t t;

  size_t n = N;
  double mass = 1.0;
  range_t range = {point_t{RMINX, RMINY, RMINZ}, point_t{RMAXX, RMAXY, RMAXZ}};
  std::cout << "Range: " << range[0] << "-" << range[1] << std::endl;

  point_t max;
  point_t min;

  for (size_t i = 0; i < n; ++i) {
    point_t p = {uniform(RMINX, RMAXX), uniform(RMINY, RMAXY),
                 uniform(RMINZ, RMAXZ)};
    auto e = t.make_entity(key_type(range, p), p, nullptr, 0, mass, 0, HMAX);
    t.insert(e);
  }

  t.cofm(t.root(), 0, false);

  ASSERT_TRUE(t.root()->mass() == n * mass);

  for (size_t i = 0; i < n; ++i) {
    auto ent = t.get(i);

    for (size_t d = 0; d < gdimension; ++d) {
      max[d] = ent->coordinates()[d] + HMAX;
      min[d] = ent->coordinates()[d] - HMAX;
    }
    auto ns = t.find_in_box(min, max, tree_geometry_t::within_box);

    set<body_holder *> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder *> s2;

    for (size_t j = 0; j < n; ++j) {
      auto ej = t.get(j);

      bool in_box = true;
      for (size_t d = 0; d < gdimension; ++d) {
        if (ej->coordinates()[d] > max[d] || ej->coordinates()[d] < min[d]) {
          in_box = false;
          break;
        }
      }

      if (in_box) {
        s2.insert(ej);
      }
    }

    ASSERT_TRUE(s1 == s2);
  }
}

#if 0
TEST(tree_topology, neighbors_box_VARIABLE) {
  tree_topology_t t;

  size_t n = N;
  double mass = 1.0;

  point_t max;
  point_t min;
  range_t range = {point_t{RMINX,RMINY,RMINZ},point_t{RMAXX,RMAXY,RMAXZ}};
  std::cout<<"Range: "<<range[0]<<"-"<<range[1]<<std::endl;

  for(size_t i = 0; i < n; ++i){
    point_t p = {uniform(RMINX, RMAXX), uniform(RMINY, RMAXY),
        uniform(RMINZ, RMAXZ)};
    auto e = t.make_entity(key_type(range,p),p,nullptr,0,mass,0,uniform(HMIN,HMAX));
    t.insert(e);
  }

  t.cofm(t.root(),0,false);

  ASSERT_TRUE(t.root()->mass() == n*mass);

  for(size_t i = 0; i < n; ++i){
    auto ent = t.get(i);

	  for(size_t d = 0; d < gdimension; ++d ){
		  max[d] = ent->coordinates()[d]+0.00001+ent->radius();
		  min[d] = ent->coordinates()[d]-0.00001-ent->radius();
	  }
    auto ns = t.find_in_box(min,max,tree_geometry_t::intersects_sphere_box);

    set<body_holder*> s1;
    s1.insert(ns.begin(), ns.end());

    set<body_holder*> s2;

    for(size_t j = 0; j < n; ++j){
      auto ej = t.get(j);

      if(tree_geometry_t::intersects_sphere_box(min,max,
        ej->coordinates(),ej->radius())){
        s2.insert(ej);
      }
	  }
    ASSERT_TRUE(s1 == s2);
  }
}
#endif
