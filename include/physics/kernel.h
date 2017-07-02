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
    // Default 1D case
    double sigma = 2./3.;
    double result = sigma/pow(h,gdimension);
    if (0.0 <= rh && rh < 1.0) {
      result *= 1.0 - (3.0/2.0) * pow(rh,2) + (3.0/4.0) * pow(rh,3); 
      return result; 
    }else if (1.0 <= rh && rh < 2.0) {
      result *= (1.0/4.0) * pow(2-rh, 3);
      return result;
    }
    return 0.0;
  } // kernel

  // Standard gradient of spline kernel
  // for 1-2-3D  
  //static
  point_t 
  cubic_spline_gradKernel(
      point_t vecP, 
      double h)
  {
    // Default 1D case
    double sigma = 2./3.;
    double coeff = sigma/pow(h,1+gdimension);
    double r = 0;
    for(int i=0;i<gdimension;++i){
      r+= vecP[i]*vecP[i];
    }
    r = sqrt(r);
    double rh = r/h;

    point_t result{};
    if (0.0 <= rh && rh < 1.0){
      result = coeff*vecP;
      result *= ((-3.0/h)+(9.0*r/(4.0*h*h)));
    }else if(1.0 <= rh && rh < 2.0){
      result = coeff*vecP;
      result *= ((-3.0/r)+(3.0/h)+(-3.0*r/(4.0*h*h)));
    }
    return result;
  } // gradKernel 

}; // kernel

#endif // _physics_kernel_h_
