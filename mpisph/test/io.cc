#include <cinchdevel.h>
#include <cinchtest.h>

#include <cmath>
#include <iostream>
#include <mpi.h>

#include "io.h"
#include "tree_colorer.h"
#include "utils.h"

using namespace std;
using namespace flecsi;
using namespace topology;

namespace flecsi {
namespace execution {
void driver(int argc, char *argv[]) {}
} // namespace execution
} // namespace flecsi

TEST(io, write_N_read) {

  srand(time(NULL));
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const char *fileprefix = "io_utest";
  const char *filename = "io_utest.h5part";

  // Generate particles, write to file, read and compare
  size_t n = 1000;
  double mass = 1.0;
  double u = 1.0;

  std::vector<body> bodies(n);

  for (size_t i = 0; i < n; ++i) {
    point_t p = {(double)rand() / (double)RAND_MAX,
                 (double)rand() / (double)RAND_MAX,
                 (double)rand() / (double)RAND_MAX};
    // Create body to check the class too
    bodies[i].set_coordinates(p);
    bodies[i].set_mass(mass);
    bodies[i].set_radius(u);
  }

  // Write that to file
  io::outputDataHDF5(bodies, fileprefix, 0, 0.);

  std::vector<body> rbodies;
  int64_t totalnbodies = 0;
  int64_t localnbodies = 0;
  // Read this file
  io::inputDataHDF5(rbodies, fileprefix, fileprefix, totalnbodies, localnbodies,
                    0);

  ASSERT_TRUE(localnbodies == n);
  ASSERT_TRUE(totalnbodies == n * size);

  // Remove the created file
  remove(filename);
}
