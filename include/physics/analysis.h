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

  void 
  display(const char * filename = NULL, bool header = false)
  {
    // TODO: output just a single line on a screen, containing iteration and time;
    //       output all scalar reductions as single clearly formatted line to a 
    //       file ("reductions.dat") as well as on the screen. When creating the 
    //       file, write a header to indicate which quantities are output in which 
    //       column (because their order and quantity may change between revisions)
    //       E.g.:
    //       -- >> example output file >> -----------------------------------------
    //       # Scalar reductions:
    //       # 1:iteration 2:time 3:energy 4:mom_x 5:mom_y 6:mom_z
    //       0  0.0   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00 
    //       10 0.1   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00 
    //       20 0.2   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00 
    //       30 0.3   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00 
    //       ...
    //       -- << end output file <<<< -------------------------------------------
    
    // For screen otuput, generate the header anyways 
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

    std::ostringstream oss_data;
    oss_data << std::setprecision(10);
    oss_data << physics::iteration <<" "<< physics::totaltime<< " ";
    for(int i = 0 ; i < gdimension ; ++i){
      oss_data<<linear_momentum[i]<<" ";
    }
    oss_data<<std::endl;

    // Output to screen 
    clog(info)<<oss_header.str();
    clog(info)<<oss_data.str();

    // If file output 
    // Header in file 
    if(filename != NULL){
      // Print header if needed, in this case erase content
      if(header == true){
        std::ofstream out(filename);
        out << oss_header.str();
        out << oss_data.str();
        out.close(); 
      }else{
        // Print the new line of data in this case append
        std::ofstream out(filename,std::ios_base::app);
        out << oss_data.str();
        out.close();
      }
    }
  }

}; // physics

#endif // _physics_analysis_h_
