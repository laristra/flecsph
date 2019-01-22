#include <cinchtest.h>

#include <mpi.h>

namespace analysis{
  enum e_conservation: size_t
  { MASS = 0 , ENERGY = 1, MOMENTUM = 2, ANG_MOMENTUM = 3 };
}
using namespace analysis;

namespace flecsi{
namespace execution{
  void mpi_init_task(const char * parameter_file);
  bool check_conservation(const std::vector<e_conservation>&);
}
}

using namespace flecsi;
using namespace execution;

TEST(sodtube, working) {
  int provided;
  MPI_Query_thread(&provided);
  ASSERT_TRUE(provided == MPI_THREAD_MULTIPLE);
  mpi_init_task("sodtube_t1_n100.par");
  ASSERT_TRUE(check_conservation({MASS,ENERGY,MOMENTUM}));
}
