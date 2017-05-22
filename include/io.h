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
 * @file io.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Simple implementation for Input Output for serial and distributed 
 */

#ifndef io_h
#define io_h 

#include <cstdlib>
#include <iostream>
#include <vector>

#include <hdf5.h>

#include "tree.h"

namespace io{
  
  void 
  inputDataTxt(
      std::vector<body*>&, 
      const char *, 
      tree_topology_t&); 

  void 
  inputDataHDF5(
      std::vector<body*>&, 
      const char *, 
      tree_topology_t&); 

  void 
  inputDataTxtRange(
      std::vector<std::pair<entity_key_t,body>>&,
      int&,
      int&,
      int,
      int,
      const char *);

  void 
  outputDataTxt();

  void 
  outputDataHDF5(
      std::vector<body*>&,
      const char*, 
      int, 
      double); 

} // namespace io

#endif // io_h

