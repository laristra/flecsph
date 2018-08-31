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

TEST(sodtube, working) {
  MPI_Init(NULL,NULL);
  mpi_init_task("sodtube_t1_n100.par");
  MPI_Finalize();
}
