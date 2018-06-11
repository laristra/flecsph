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
 * @file physics.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic physics implementation
 */

#ifndef _physics_kernel_h_
#define _physics_kernel_h_

#include <vector>

#include "tree.h"

namespace kernel{

  // Standard spline kernel
  // for 1-2-3D
  //static 
  double 
  cubic_spline_kernel(
      double r, 
      double h)
  {
    double rh = r/h;
    // Default 3D case
    double sigma = 1.0/(M_PI*h*h*h);
    //if(gdimension == 2){
    //  sigma = 10./(7.*M_PI*h*h);
    //}
    double result = 0.; 
    if (0.0 <= rh && rh < 1.0) {
      result = 1.0 - 1.5 * rh*rh*(1-.5*rh);
      result *= sigma;  
    }else if (1.0 <= rh && rh < 2.0) {
      result *= 0.25 * (2-rh)*(2-rh)*(2-rh);
      result *= sigma; 
    }
    return result;
  } // kernel

  // Standard gradient of spline kernel
  // for 1-2-3D  
  //static
  point_t 
  cubic_spline_gradKernel(
      point_t vecP, 
      double h)
  {
    // Default 3D case
    double sigma = 1.0/(M_PI*h*h*h);
    //if(gdimension == 2){
    //  sigma = 10./(7.*M_PI*h*h);
    //}
    //double coeff = sigma; ///pow(h,1+gdimension);
    
    // Distance 
    double r = 0;
    for(size_t i=0;i<gdimension;++i){
      r += vecP[i]*vecP[i];
    }
    r = sqrt(r);
    double rh = r/h;

    point_t result{};
    if (0.0 <= rh && rh < 1.0){
      result = sigma*vecP;
      result *= -3.0 * rh * (1-.75*rh)/h/r;
    }else if(1.0 <= rh && rh < 2.0){
      result = sigma*vecP;
      result *= -.75*(2-rh)*(2-rh)/h/r;
    }
    return result;
  } // gradKernel

};

#endif // _physics_kernel_h_
