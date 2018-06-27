#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>

#include <mpi.h>

namespace flecsi{
namespace execution{
  void mpi_init_task(int numberiterations);

}
}

using namespace flecsi;
using namespace execution;

TEST(sedov, working) {
  MPI_Init(NULL,NULL);
  mpi_init_task(25);
  MPI_Finalize();
}
