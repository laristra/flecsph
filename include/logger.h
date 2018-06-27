/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

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
 * @file logger.h
 * @author Oleg Korobkin
 * @date June 2018
 * @brief Logger class for output
 *        Inspired by: http://www.drdobbs.com/cpp/logging-in-c/201804215
 */

#ifndef logger_h
#define logger_h
#define LOGGER (*logger::out)

namespace logger {
  std::ostream *original_cin;
  std::ofstream ofs_null;
  std::ostream *out;

  void init() {
    int my_rank; 
    ofs_null.open ("/dev/null");
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    original_cin = std::cin.tie ();

    if (my_rank > 0)
      // redirect stdout to /dev/null
      original_cin = std::cin.tie (&ofs_null);

    out = std::cin.tie();

  }
  void finalize() {
    std::cin.tie (original_cin);
    ofs_null.close();
  }
} //namespace logger

#endif // logger_h

