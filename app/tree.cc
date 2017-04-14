#include "flecsi.h"
#include "flecsi/execution/execution.h"
#include "flecsi/concurrency/thread_pool.h"

#if FLECSI_RUNTIME_MODEL == FLECSI_RUNTIME_MODEL_mpilegion
  #include <mpi.h>
  #include <legion.h>
#endif

int main(int argc, char * argv[]){

  #if FLECSI_RUNTIME_MODEL == FLECSI_RUNTIME_MODEL_mpilegion
    MPI_Init(&argc,&argv);
  #endif

  auto retval = flecsi::execution::context_t::instance().initialize(argc,argv);
  return retval;
  
  #if FLECSI_RUNTIME_MODEL == FLECSI_RUNTIME_MODEL_mpilegion
    MPI_Finalize();
  #endif

}

