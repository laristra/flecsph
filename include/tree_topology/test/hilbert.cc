#include <cinchdevel.h>
#include <cinchtest.h>

#include <cmath>
#include <iostream>
#include <mpi.h>

#include "../hilbert_id.h"

using namespace std;
using namespace flecsi;
using namespace topology;

namespace flecsi {
namespace execution {
void driver(int argc, char *argv[]) {}
} // namespace execution
} // namespace flecsi

const double tol = 1.0e-6;

template <typename T, size_t D>
bool operator==(const point__<T, D> &a, const point__<T, D> &b) {
  for (size_t d = 0; d < D; ++d) {
    if (a[d] != b[d])
      return false;
  }
  return true;
}

using point_2d = point__<double, 2>;
using hilbert_2d = hilbert_id<uint64_t, 2>;
const std::array<point_2d, 2> range_2d{point_2d{0., 0.}, point_2d{1., 1.}};

TEST(hilbert_2d, coordinates) {
  size_t dimension = 2;

  for (int t = 0; t < 100; ++t) {
    point__<double, 2> pt(
        (double)rand() / (double)RAND_MAX * (range_2d[1][0] - range_2d[0][0]) +
            range_2d[0][0],
        (double)rand() / (double)RAND_MAX * (range_2d[1][1] - range_2d[0][1]) +
            range_2d[0][1]);
    hilbert_2d ent = hilbert_2d(range_2d, pt);
    point__<double, 2> ptent;
    ent.coordinates(range_2d, ptent);
    std::cout << pt << "==" << ptent << std::endl;
    ASSERT_TRUE(abs(pt[0] - ptent[0]) < tol);
    ASSERT_TRUE(abs(pt[1] - ptent[1]) < tol);
  }
}

using point_3d = point__<double, 3>;
using hilbert_3d = hilbert_id<uint64_t, 3>;
const std::array<point_3d, 2> range_3d{point_3d{0., 0., 0.},
                                       point_3d{1., 1., 1.}};

TEST(hilbert_3d, range_branch) {
  size_t dimension = 3;

  for (int t = 0; t < 2; ++t) {
    point__<double, 3> pt(
        (double)rand() / (double)RAND_MAX * (range_3d[1][0] - range_3d[0][0]) +
            range_3d[0][0],
        (double)rand() / (double)RAND_MAX * (range_3d[1][1] - range_3d[0][1]) +
            range_3d[0][1],
        (double)rand() / (double)RAND_MAX * (range_3d[1][2] - range_3d[0][2]) +
            range_3d[0][2]);
    hilbert_3d ent = hilbert_3d(range_3d, pt);
    point__<double, 3> ptent;
    ent.coordinates(range_3d, ptent);
    std::cout << pt << "==" << ptent << std::endl;
    ASSERT_TRUE(abs(pt[0] - ptent[0]) < tol);
    ASSERT_TRUE(abs(pt[1] - ptent[1]) < tol);
    ASSERT_TRUE(abs(pt[2] - ptent[2]) < tol);
  }
}
