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

#include "physics.h"
#include "logger.h"

namespace analysis{

  point_t linear_momentum;
  double JULIEN_ENERGY; 

  double lambda_julien = 4.5;

  // Calculate simple linear momentum for checking momentum conservation
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
  compute_JULIEN_ENERGY(
    std::vector<body_holder*>& bodies)
  {
    JULIEN_ENERGY = 0;
    for(auto nbh: bodies){
      JULIEN_ENERGY += lambda_julien* 
      nbh->getBody()->getDensity() / nbh->getBody()->getMass()
      * nbh->getBody()->getInternalenergy() * nbh->getBody()->getSmoothinglength(); 
    }
  }

  void 
  display()
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
    //
    LOGGER << std::endl << "Analysis: " <<std::endl
           << "Linear momentum: " << linear_momentum << std::endl
           << "JULIEN_ENERGY: " << JULIEN_ENERGY << std::endl<<std::endl;
  }

}; // physics

#endif // _physics_analysis_h_
