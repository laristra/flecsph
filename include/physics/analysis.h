/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@ @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@///// /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@      /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@ /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////  /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@      /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@      /@@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file analysis.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic analysis
 */

#ifndef _physics_analysis_h_
#define _physics_analysis_h_

#include <vector>

//#include "physics.h"

namespace analysis{

  point_t linear_momentum;

  /**
   * @brief      Compute the linear momentum,
   *
   * @param      bodies  Vector of all the local bodies
   * @param      total   Total linear momentum
   */
  void
  compute_lin_momentum(
      std::vector<body_holder*>& bodies)
  {
    linear_momentum = {0};
    for(auto nbh: bodies) {
      linear_momentum += nbh->getBody()->getLinMomentum();
    }
  }


  /**
   * @brief Rolling screen output
   */
  void
  screen_output() 
  {
    const int screen_length = 40;  
    (physics::iteration-1)%screen_length  ||
    clog_one(info)<< "#-- iteration:               time:" <<std::endl;
    clog_one(info)
      << std::setw(14) << physics::iteration
      << std::setw(20) << std::scientific << std::setprecision(12)
      << physics::totaltime << std::endl;
  }


  /**
   * @brief Periodic file output
   *
   * Outputs all scalar reductions as single formatted line to a file. When
   * creating the file, writes a header to indicate which quantities are
   * output in which column (because their order and quantity may change
   * between revisions)
   * E.g.:
   * -- >> example output file >> -----------------------------------------
   * # Scalar reductions:
   * # 1:iteration 2:time 3:energy 4:mom_x 5:mom_y 6:mom_z
   * 0  0.0   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * 10 0.1   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * 20 0.2   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * 30 0.3   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * ...
   * -- << end output file <<<< -------------------------------------------
   */
  void
  scalar_output(const char * filename = NULL)
  {
    static bool first_time = true;

    if (first_time) {
      // Generate and output the header
      std::ostringstream oss_header;
      oss_header
        << "# Scalar reductions: " <<std::endl
        << "# 1:iteration 2:time 3:mom_x ";

      // The momentum depends on dimension
      if(gdimension > 1){
        oss_header<<"4:mom_y ";
      }
      if(gdimension == 3){
        oss_header<<"5:mom_z ";
      }
      oss_header<<std::endl;

      std::ofstream out(filename);
      out << oss_header.str();
      out.close();

      first_time = false;
    }

    std::ostringstream oss_data;
    oss_data << std::setw(14) << physics::iteration
      << std::setw(20) << std::scientific << std::setprecision(12)
      << physics::totaltime;
    for(int i = 0 ; i < gdimension ; ++i){
      oss_data
        << std::setw(20) << std::scientific << std::setprecision(12)
        << linear_momentum[i] <<" ";
    }
    oss_data << std::endl;

    // Open file in append mode
    std::ofstream out(filename,std::ios_base::app);
    out << oss_data.str();
    out.close();

  } // scalar output

}; // physics

#endif // _physics_analysis_h_
