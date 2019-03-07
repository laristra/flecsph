#include <cinchdevel.h>
#include <cinchtest.h>

#include <cmath>
#include <iostream>
#include <mpi.h>

#include "../filling_curve.h"

using namespace std;
using namespace flecsi;

namespace flecsi {
namespace execution {
void driver(int argc, char *argv[]) {}
} // namespace execution
} // namespace flecsi

template <typename T, size_t D>
bool operator==(const point__<T, D> &a, const point__<T, D> &b) {
  for (size_t d = 0; d < D; ++d) {
    if (a[d] != b[d])
      return false;
  }
  return true;
}

using point_2d = point__<double, 2>;
using morton_2d = morton_curve<2, uint64_t>;
const std::array<point_2d, 2> range_2d{point_2d{0., 0.}, point_2d{1., 1.}};

TEST(morton_2d, range_branch) {
  size_t dimension = 2;

  morton_2d root = morton_2d::root();
  auto root_range = root.range(range_2d);
  ASSERT_TRUE(root_range[0] == range_2d[0] && root_range[1] == range_2d[1]);
  // std::cout<<"rt  "<<"k="<<root<<" ["<<root_range[0]<<";"<<
  //  root_range[1]<<"]"<<std::endl;

  // Check the first children for this dimension
  for (size_t d = 0; d < 1 << dimension; ++d) {
    // Create child key
    morton_2d child = root;
    child.push(d);
    // Compute its range
    auto rg = child.range(range_2d);
    // std::cout<<d<<" "<<"k="<<child<<" ["<<rg[0]<<";"<<rg[1]<<"]"<<std::endl;
    // entity test, should be inside the range with its key
    point__<double, 2> pt(
        (double)rand() / (double)RAND_MAX * (rg[1][0] - rg[0][0]) + rg[0][0],
        (double)rand() / (double)RAND_MAX * (rg[1][1] - rg[0][1]) + rg[0][1]);
    morton_2d ent = morton_2d(range_2d, pt);
    // std::cout<<"Entity key: "<< pt << " k="<<ent<<std::endl;
    // Assert that the range and the generated entity correspond to this
    // subspace
    ent.truncate(1);
    ASSERT_TRUE(ent == child);
  }

  root.push(2);
  // Check the first children for this dimension
  for (size_t d = 0; d < 1 << dimension; ++d) {
    // Create child key
    morton_2d child = root;
    child.push(d);
    // Compute its range
    auto rg = child.range(range_2d);
    // std::cout<<d<<" "<<"k="<<child<<" ["<<rg[0]<<";"<<rg[1]<<"]"<<std::endl;
    // entity test, should be inside the range with its key
    for (int t = 0; t < 50; ++t) {
      point__<double, 2> pt(
          (double)rand() / (double)RAND_MAX * (rg[1][0] - rg[0][0]) + rg[0][0],
          (double)rand() / (double)RAND_MAX * (rg[1][1] - rg[0][1]) + rg[0][1]);
      morton_2d ent = morton_2d(range_2d, pt);
      // std::cout<<"Entity key: "<< pt << " k="<<ent<<std::endl;
      // Assert that the range and the generated entity correspond to this
      // subspace
      ent.truncate(2);
      ASSERT_TRUE(ent == child);
    }
  }
}

using point_3d = point__<double, 3>;
using morton_3d = morton_curve<3, uint64_t>;
const std::array<point_3d, 2> range_3d{point_3d{0., 0., 0.},
                                       point_3d{1., 1., 1.}};

TEST(morton_3d, range_branch) {
  size_t dimension = 3;

  morton_3d root = morton_3d::root();
  auto root_range = root.range(range_3d);
  ASSERT_TRUE(root_range[0] == range_3d[0] && root_range[1] == range_3d[1]);
  // std::cout<<"rt  "<<"k="<<root<<" ["<<root_range[0]<<";"<<
  //  root_range[1]<<"]"<<std::endl;

  // Check the first children for this dimension
  for (size_t d = 0; d < 1 << dimension; ++d) {
    // Create child key
    morton_3d child = root;
    child.push(d);
    // Compute its range
    auto rg = child.range(range_3d);
    // std::cout<<d<<" "<<"k="<<child<<" ["<<rg[0]<<";"<<rg[1]<<"]"<<std::endl;
    // entity test, should be inside the range with its key
    point__<double, 3> pt(
        (double)rand() / (double)RAND_MAX * (rg[1][0] - rg[0][0]) + rg[0][0],
        (double)rand() / (double)RAND_MAX * (rg[1][1] - rg[0][1]) + rg[0][1],
        (double)rand() / (double)RAND_MAX * (rg[1][2] - rg[0][2]) + rg[0][2]);
    morton_3d ent = morton_3d(range_3d, pt);
    // std::cout<<"Entity key: "<< pt << " k="<<ent<<std::endl;
    // Assert that the range and the generated entity correspond to this
    // subspace
    ent.truncate(1);
    ASSERT_TRUE(ent == child);
  }

  root.push(2);
  // Check the first children for this dimension
  for (size_t d = 0; d < 1 << dimension; ++d) {
    // Create child key
    morton_3d child = root;
    child.push(d);
    // Compute its range
    auto rg = child.range(range_3d);
    // std::cout<<d<<" "<<"k="<<child<<" ["<<rg[0]<<";"<<rg[1]<<"]"<<std::endl;
    // entity test, should be inside the range with its key
    for (int t = 0; t < 50; ++t) {
      point__<double, 3> pt(
          (double)rand() / (double)RAND_MAX * (rg[1][0] - rg[0][0]) + rg[0][0],
          (double)rand() / (double)RAND_MAX * (rg[1][1] - rg[0][1]) + rg[0][1],
          (double)rand() / (double)RAND_MAX * (rg[1][2] - rg[0][2]) + rg[0][2]);
      morton_3d ent = morton_3d(range_3d, pt);
      // std::cout<<"Entity key: "<< pt << " k="<<ent<<std::endl;
      // Assert that the range and the generated entity correspond to this
      // subspace
      ent.truncate(2);
      ASSERT_TRUE(ent == child);
    }
  }
}
