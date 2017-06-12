/*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file main_driver.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Specialization and Main driver used in FleCSI. 
 * The Specialization Driver is normally used to register data and the main 
 * code is in the Driver.  
 */

#include <iostream>
#include <numeric> // For accumulate
#include <iostream>

#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

#include "flecsi/io/simple_definition.h"
#include "flecsi/coloring/tree_colorer.h"
//#include "flecsi/coloring/parmetis_colorer.h"
//#include "flecsi/coloring/mpi_communicator.h"
//#include "flecsi/topology/closure_utils.h"
#include "flecsi/utils/set_utils.h"

#include <bodies_system.h>
#include "io.h"
#include "physics.h"

namespace flecsi{
namespace execution{

void
mpi_init_task(
    int startiteration)
{
  using entry_info_t = flecsi::coloring::entity_info_t;

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int64_t localnbodies, totalnbodies; 
  std::vector<std::pair<entity_key_t,body>> localbodies; 

  // Read the particles 
  io::inputDataHDF5(localbodies,"hdf5_sodtube.h5part",
      totalnbodies,localnbodies);

  // Compute the keys 

  // Compute the ghosts, shared, exclusive

}

flecsi_register_task(mpi_init_task,mpi,index);

void 
specialization_driver(
    int argc, 
    char * argv[])
{
  // Default start at iteration 0
  int startiteration = 0;
  if(argc == 2){
    startiteration = atoi(argv[1]);
  }
  std::cout << "In user specialization_driver" << std::endl;
  
  flecsi_execute_task(mpi_init_task,mpi,index,startiteration); 
} // specialization driver

void 
driver(
    int argc,  
    char * argv[])
{
  std::cout << "In user driver" << std::endl;

  // Do the first loop over the data 
} // driver


} // namespace
} // namespace


