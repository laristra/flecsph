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

TEST(noh, working) {
  //MPI_Init(NULL,NULL);
  mpi_init_task("noh_nx20.par");
  //MPI_Finalize();
}
