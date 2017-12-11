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


  double kernel_2D_fac = 7/(4*M_PI);
  double grad_kernel_2D_fac = -35./(4.*M_PI);

  double kernel_3D_fac = 21/(16*M_PI);
  double grad_kernel_3D_fac = -147./(16.*M_PI);

  
  double quintic_wendland_2D(
      double r, 
      double h)
  {
    double q = r/h;
    return kernel_2D_fac/(h*h)*pow(1.-0.5*q,4)*(2.*q+1.); 
  }

  point_t quintic_wendland_2D_gradient(
      point_t vec, 
      double h)
  {
    double r = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
    double q = r/h;
    return grad_kernel_2D_fac/(h*h*h) * q * pow(1.-0.5*q,3)/r * vec;
  } 

  double quintic_wendland_3D(
      double r, 
      double h)
  {
    double q = r/h;
    return kernel_3D_fac/(h*h*h)*pow(1.-0.5*q,4)*(2.*q+1.); 
  }

  point_t quintic_wendland_3D_gradient(
      point_t vec, 
      double h)
  {
    double r = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
    double q = r/h;
    return grad_kernel_3D_fac/(h*h*h*h) * q * pow(1.-0.5*q,3)/r * vec;
  } 


}; // kernel

#endif // _physics_kernel_h_
