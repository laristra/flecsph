#include <cinchdevel.h>
#include <cinchtest.h>

#include <cmath>
#include <iostream>
#include <mpi.h>

#include "bodies_system.h"

using namespace std;
using namespace flecsi;
using namespace topology;

namespace flecsi {
namespace execution {
void driver(int argc, char *argv[]) {}
} // namespace execution
} // namespace flecsi

TEST(body_system, write_range_read) {

  const char *fileprefix = "io_test";

  body_system<double, gdimension> bs;
  bs.read_bodies(fileprefix, fileprefix, 0);

  double h = bs.getSmoothinglength();
  ASSERT_TRUE(h == 0.05);
  std::array<point_t, 2> range = bs.getRange();
  for (size_t i = 0; i < gdimension; ++i) {
    ASSERT_TRUE(range[0][i] == -0.05);
    ASSERT_TRUE(fabs(range[1][i] - 0.50) < 1.0e-15);
  }

  bs.write_bodies(fileprefix, 0, 0);
}
