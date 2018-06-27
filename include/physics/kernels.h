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

namespace kernels{

  /**
   * @brief      Compute the size of the vector 
   * needed in the computation of the gradients
   *
   * @param      p     The vector, in point_t shape
   *
   * @return     The size of the vector
   */
  double vector_size(
    point_t& p){
    double result = 0.;
    for(int i = 0 ; i < gdimension; ++i){
      result += p[i]*p[i];
    }
    return sqrt(result);
  }

  // Coefficient for the kernels in 1d, 2d and 3d without h 
  // h depends of the dimension and is added in the kernel
  const double cubic_spline_sigma[3] = {2./3.,10./(7.*M_PI),1./M_PI};
  const double gaussian_sigma[3] = {1./pow(M_PI,.5),1./M_PI,1./pow(M_PI,1.5)};
  const double quintic_spline_sigma[3] = {1./120.,7./(478.*M_PI),3./(359.*M_PI)};
  const double wendland_quintic_sigma[3] = {5./8.,7./(4.*M_PI),21./(16.*M_PI)};

/*============================================================================*/
/*   Cubic spline                                                             */
/*============================================================================*/
  /**
   * @brief      Cubic splein kernel
   * From Monaghan/92
   *
   * @param[in]  r     Distance between the particles 
   * @param[in]  h     Smoothing length 
   *
   * @return     Contribution from the particle 
   */
  double 
  cubic_spline(
      double r, 
      double h)
  {
    double rh = r/h;
    // \TODO need to use solution based on template 
    double sigma = cubic_spline_sigma[gdimension-1]/pow(h,gdimension);
    double result = 0.;

    if (0.0 <= rh && rh <= 1.0) {
      result = 1.0 - 1.5*rh*rh + .75*rh*rh*rh;
      result *= sigma;  
    }else if (1.0 < rh && rh <= 2.0) {
      result = 0.25 * (2-rh)*(2-rh)*(2-rh);
      result *= sigma; 
    }
    return result;
  } // kernel

  /**
   * @brief      Gradient of cubic spline kernel
   * From Monaghan/92
   *
   * @param[in]  vecP  The vector pab = pa - pb 
   * @param[in]  h     The smoothing length 
   *
   * @return     Contribution from the particle 
   */
  point_t 
  gradient_cubic_spline(
      point_t vecP, 
      double h)
  {
    double sigma = cubic_spline_sigma[gdimension-1]/pow(h,gdimension+1);
    // Compute distance of particles 
    double r = vector_size(vecP);
    // Normalize vector 
    point_t eab = vecP / r;
    double rh = r/h;

    point_t result{};
    if (0.0 <= rh && rh <= 1.0){
      result = sigma*eab;
      result *= -3.0*rh + 9./4.*rh*rh;
    }else if(1.0 < rh && rh <= 2.0){
      result = sigma*eab;
      result *= -.75*(2-rh)*(2-rh);
    }
    return result;
  } // gradKernel 

/*============================================================================*/
/*   Gaussian                                                                 */
/*============================================================================*/
  /**
   * @brief      Gaussian kernel 
   * From Liu/2010
   *
   * @param[in]  r     Distance between the particles 
   * @param[in]  h     Smoothing length 
   *
   * @return     Contribution from the particle 
   */
  double 
  gaussian(
    double r, 
    double h)
  {
    double rh = r/h;
    double sigma = gaussian_sigma[gdimension-1]/pow(h,gdimension);
    double result = 0.;
    if(rh <= 3.){
      result = sigma*exp(-rh*rh);
    }
    return result; 
  }

  /**
   * @brief      Gradient of gaussian kernel
   * From Liu/2010
   *
   * @param[in]  vecP  The vector pab = pa - pb 
   * @param[in]  h     The smoothing length 
   *
   * @return     Contribution from the particle 
   */
  point_t 
  gradient_gaussian(
    point_t vecP,
    double h)
  {
    double sigma = gaussian_sigma[gdimension-1]/pow(h,gdimension+1);
    // Compute distance of particles 
    double r = vector_size(vecP);
    // Normalize vector 
    point_t eab = vecP / r;
    double rh = r/h;

    point_t result{};
    if (rh <= 3.){
      result = sigma*eab; 
      result *= -2.*rh*exp(-rh*rh);
    }
    return result;
  }

/*============================================================================*/
/*   Quintic spline                                                           */
/*============================================================================*/
    /**
   * @brief      Quintic spline kernel 
   * From Liu/2010
   *
   * @param[in]  r     Distance between the particles 
   * @param[in]  h     Smoothing length 
   *
   * @return     Contribution from the particle 
   */
  double 
  quintic_spline(
    double r, 
    double h)
  {
    double rh = r/h;
    double sigma = quintic_spline_sigma[gdimension-1]/pow(h,gdimension);
    double result = 0.;
    if(0 <= rh && rh <= 1){
      result = sigma*(pow(3-rh,5)-6*pow(2-rh,5)+15*pow(1-rh,5));
    }else if(1 < rh && rh <= 2){
      result = sigma*(pow(3-rh,5)-6*pow(2-rh,5));
    }else if(2 < rh && rh <= 3){
      result = sigma*(pow(3-rh,5));
    }
    return result; 
  }

  /**
   * @brief      Gradient of quintic spline kernel
   * From Liu/2010
   *
   * @param[in]  vecP  The vector pab = pa - pb 
   * @param[in]  h     The smoothing length 
   *
   * @return     Contribution from the particle 
   */
  point_t 
  gradient_quintic_spline(
    point_t vecP,
    double h)
  {
    double sigma = quintic_spline_sigma[gdimension-1]/pow(h,gdimension+1);
    // Compute distance of particles 
    double r = vector_size(vecP);
    // Normalize vector 
    point_t eab = vecP / r;
    double rh = r/h;

    point_t result{};
    if(0 <= rh && rh <= 1){
      result = sigma*eab;
      result *= (-5.*pow(3-rh,4)+30.*pow(2-rh,4)-75.*pow(1-rh,4));
    }else if(1 < rh && rh <= 2){
      result = sigma*eab;
      result *= (-5.*pow(3-rh,4)+30.*pow(2-rh,4));
    }else if(2 < rh && rh <= 3){
      result = sigma*eab;
      result *= (-5.*(pow(3-rh,4)));
    }
    return result; 
  }

/*============================================================================*/
/*   Wendland quintic                                                         */
/*============================================================================*/
    /**
   * @brief      Wendland quintic
   * \TODO add the ref
   *
   * @param[in]  r     Distance between the particles 
   * @param[in]  h     Smoothing length 
   *
   * @return     Contribution from the particle 
   */
  double 
  wendland_quintic(
    double r, 
    double h)
  {
    double rh = r/h;
    double sigma = wendland_quintic_sigma[gdimension-1]/pow(h,gdimension);
    double result = 0.;
    // Different cases for 1D and 2D/3D 
    if(gdimension == 1){
      if(0 <= rh && rh <= 2){
        result = sigma*(pow(1-rh/2.,3)*(1.5*rh+1));
      }
    }else{
      if(0 <= rh && rh <= 2){
        result = sigma*(pow(1-rh/2.,4)*(2 *rh+1));
      }
    }
    return result; 
  }

  /**
   * @brief      Gradient of quintic spline kernel
   * \TODO add the ref
   *
   * @param[in]  vecP  The vector pab = pa - pb 
   * @param[in]  h     The smoothing length 
   *
   * @return     Contribution from the particle 
   */
  point_t 
  gradient_wendland_quintic(
    point_t vecP,
    double h)
  {
    double sigma = wendland_quintic_sigma[gdimension-1]/pow(h,gdimension+1);
    // Compute distance of particles 
    double r = vector_size(vecP);
    // Normalize vector 
    point_t eab = vecP / r;
    double rh = r/h;

    point_t result{};
    if(gdimension == 1){
      if(0 <= rh && rh <= 2){
        result = sigma*eab; 
        result *= -3.*rh*pow(1.-rh/2.,2);
      }
    }else{
      if(0 <= rh && rh <= 2){
        result = sigma*eab; 
        result *= -5.*rh*pow(1.-rh/2.,3);
      }
    }
    return result; 
  }


}; // kernel

#endif // _physics_kernel_h_
