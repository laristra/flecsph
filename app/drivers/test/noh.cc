#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>

#include <mpi.h>

namespace flecsi{
namespace execution{
  void mpi_init_task(const char * parameter_file);
}
}

using namespace flecsi;
using namespace execution;

TEST(sedov, working) {
  mpi_init_task("noh_sqn100.par");
}
