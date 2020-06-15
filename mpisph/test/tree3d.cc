#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <log.h>

#include "default_physics.h"
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

std::ostream &
operator<<(std::ostream & ostr, const key_type id) {
  id.output_(ostr);
  return ostr;
}

double
uniform() {
  return double(rand()) / RAND_MAX;
}

double
uniform(double a, double b) {
  return a + (b - a) * uniform();
}

namespace flecsi {
namespace execution {
void
driver(int argc, char * argv[]) {}
} // namespace execution
} // namespace flecsi

TEST(tree_topology, neighbors_sphere_NORMAL) {
  MPI_Init(nullptr, nullptr);
  tree_topology_t t;

  size_t n = N;
  double mass = 1.0;
  range_t range = {point_t{RMINX, RMINY, RMINZ}, point_t{RMAXX, RMAXY, RMAXZ}};
  std::cout << "Range: " << range[0] << "-" << range[1] << std::endl;

  t.set_range(range);

  std::vector<body> entities;

  for(size_t i = 0; i < n; ++i) {
    point_t p = {
      uniform(RMINX, RMAXX), uniform(RMINY, RMAXY), uniform(RMINZ, RMAXZ)};
    t.entities().push_back(body{});
    t.entities().back().set_coordinates(p);
    t.entities().back().set_mass(mass);
    t.entities().back().set_radius(HMAX);
  }

  t.compute_keys();

  std::sort(
    t.entities().begin(), t.entities().end(), [](auto & left, auto & right) {
      if(left.key() < right.key()) {
        return true;
      }
      if(left.key() == right.key()) {
        return left.id() < right.id();
      }
      return false;
    }); // sort

  t.build_tree(physics::compute_cofm);

  ASSERT_TRUE(t.get_node(t.root())->mass() == n * mass);

  for(size_t i = 0; i < n; ++i) {
    auto ent = &(t.entities()[i]);

    // std::cout<<"Entity "<<i+1<<"/"<<n<<" = "<<ent->key();
    auto ns =
      t.find_in_radius(ent->coordinates(), HMAX, tree_geometry_t::within);
    // std:cout<<" -> "<<ns.size()<<" nbrs.done"<<std::endl;

    set<body *> s1;
    s1.insert(ns.begin(), ns.end());

    set<body *> s2;

    for(size_t j = 0; j < n; ++j) {
      auto ej = &(t.entities()[j]);

      if(distance(ent->coordinates(), ej->coordinates()) < HMAX) {
        s2.insert(ej);
      }
    }

    ASSERT_TRUE(s1 == s2);
  }
  MPI_Finalize();
}

#if 0
TEST(tree_topology, neighbors_sphere_VARIABLE) {
  tree_topology_t t;

  size_t n = N;
  double mass = 1.0;
  range_t range = {point_t{RMINX, RMINY, RMINZ}, point_t{RMAXX, RMAXY, RMAXZ}};
  std::cout << "Range: " << range[0] << "-" << range[1] << std::endl;
  //key_type::set_range(range);

  for (size_t i = 0; i < n; ++i) {
    point_t p = {uniform(RMINX, RMAXX), uniform(RMINY, RMAXY),
                 uniform(RMINZ, RMAXZ)};
    auto e = t.make_entity(key_type(range,p), p, nullptr, 0, mass, 0,
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
      if (dist <= std::max(ent->radius(),ej->radius())) {
        s2.insert(ej);
      }
    }

    std::cout<<s1.size()<<" "<<s2.size()<<std::endl;
    ASSERT_TRUE(s1 == s2);
  }
}
#endif
