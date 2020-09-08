#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <log.h>
#include <mpi.h>

#include "../../tree.h"
#include "default_physics.h"

using namespace flecsi;
using namespace topology;

namespace flecsi {
namespace execution {
void
driver(int, char **) {}
} // namespace execution
} // namespace flecsi

const size_t dimension = gdimension;
// using range_t = std::array<point_t,2>;

TEST(tree, add_entities) {
  MPI_Init(nullptr, nullptr);
  range_t range{point_t(0., 0., 0.), point_t(1., 1., 1.)};

  tree_topology_t * tree;
  tree = new tree_topology_t();
  tree->set_range(range);

  size_t nbodies = 1000;
  size_t id = 0;
  for(size_t i = 0; i < nbodies; ++i) {
    tree->entities().push_back(body{});
    tree->entities()[i].set_coordinates(
      point_t((double)rand() / (double)RAND_MAX,
        (double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX));
    tree->entities()[i].set_mass((double)rand() / (double)RAND_MAX);
    tree->entities()[i].set_id(id++);
    tree->entities()[i].set_radius((double)rand() / (double)RAND_MAX);
  }

  tree->compute_keys();

  std::sort(tree->entities().begin(), tree->entities().end(),
    [](auto & left, auto & right) {
      if(left.key() < right.key()) {
        return true;
      }
      if(left.key() == right.key()) {
        return left.id() < right.id();
      }
      return false;
    }); // sort

  tree->build_tree(physics::compute_cofm);

  // Destroy the tree
  delete tree;
  MPI_Finalize();
}
