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
    std::cout<<std::endl<<"Analysis: "<<std::endl;
    std::cout<<"Linear momentum: "<<linear_momentum<<std::endl;
    std::cout<<"JULIEN_ENERGY: "<<JULIEN_ENERGY<<std::endl<<std::endl;
  }

}; // physics

#endif // _physics_analysis_h_
