#ifndef tree_driver_h
#define tree_driver_h

#include <iostream>
#include <numeric> // For accumulate
#include <iostream>

#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

#include "physics.h"
#include "io.h"

namespace flecsi{
namespace execution{

void
mpi_task(){
  std::cout << "Hi from task !" << std::endl;
}

flecsi_register_task(mpi_task,mpi,index);

void 
specialization_driver(int argc, char * argv[]){
  std::cout << "In user specialization_driver" << std::endl;
  flecsi_execute_task(mpi_task,mpi,index); 
} // specialization driver

void 
driver(int argc, char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


#endif // tree_driver_h
