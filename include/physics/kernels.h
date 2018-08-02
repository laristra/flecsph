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
  double vector_norm(
    const point_t& p){
    if (gdimension == 3)
      return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    else if (gdimension == 2)
      return sqrt(p[0]*p[0] + p[1]*p[1]);
    else
      return abs(p[0]);
  }

  // Coefficient for the kernels in 1d, 2d and 3d without h 
  // h depends of the dimension and is added in the kernel
  const double cubic_spline_sigma[3] = {4./3.,40./(7.*M_PI),8./M_PI};
  const double gaussian_sigma[3] = {
      1.69260614115414981387661988500545775405505599758743902950990362363536L,
      2.86514256233641438778811727710055595794190432955152715087643977681460L,
      4.85098600188835377710867224691152783462583160522829280247414806262620L};
  const double quintic_spline_sigma[3] = {1./40.,63./(478.*M_PI),81./(359.*M_PI)};
  const double wendland_quintic_sigma[3] = {1.25,7./M_PI,21./(2.*M_PI)};

/*============================================================================*/
/*   Cubic spline                                                             */
/*============================================================================*/
  /**
   * @brief      Cubic spline kernel
   * From Monaghan/92
   *
   * @param[in]  r     Distance between the particles 
   * @param[in]  h     Smoothing length 
   *
   * @return     Contribution from the particle 
   */
  double 
  cubic_spline(
      const double r, 
      const double h)
  {
    double rh = 2.*r/h;

    // \TODO need to use solution based on template 
    double result = 0.;

    if (rh < 2.0) {
      if (rh < 1.0) 
        result = 1.0 - 1.5*rh*rh + .75*rh*rh*rh;
      else
        result = 0.25*(2 - rh)*(2 - rh)*(2 - rh);

      result *= cubic_spline_sigma[gdimension-1]
              / pow(h,gdimension);
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
      const point_t & vecP, 
      const double h)
  {
    double r = vector_norm(vecP);
    double rh = 2.*r/h;

    point_t result = 0.0;
    if (rh < 2.0) {
      double dWdr;
      double sigma = 2.*cubic_spline_sigma[gdimension-1]
                   / pow(h,gdimension+1);
      if (rh < 1.0)
        dWdr = -3.0*rh + 9./4.*rh*rh;
      else
        dWdr = -.75*(2-rh)*(2-rh);

      result = vecP*sigma*dWdr/r;
    }
    return result;

  } // gradKernel 

/*============================================================================*/
/*   Gaussian                                                                 */
/*============================================================================*/
  /**
   * @brief      Gaussian kernel: W(r,h) = exp( -[3r/h]^2 ) 
   * From Liu/2010
   *
   * @param[in]  r     Distance between the particles 
   * @param[in]  h     Smoothing length 
   *
   * @return     Contribution from the particle 
   */
  double 
  gaussian(
    const double r, 
    const double h)
  {
    double rh = 3.*r/h;

    double result = 0.;
    if(rh <= 3.){
      double sigma = gaussian_sigma[gdimension-1]
                   / pow(h,gdimension);
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
    const point_t & vecP,
    const double h)
  {
    double r = vector_norm(vecP);
    double rh = 3.*r/h;

    point_t result = 0.0;
    if (rh < 3.0) {
      double sigma = 3.*gaussian_sigma[gdimension-1]
                   / pow(h,gdimension+1);
      double dWdr = -2.*rh*exp(-rh*rh);
      result = vecP*sigma*dWdr/r;
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
    const double r, 
    const double h)
  {
    double rh = 3.*r/h;

    double result = 0.;
    if(rh < 3.) {
      result = pow(3-rh,5);
      if (rh < 2.)
        result += -6*pow(2-rh,5);

      if (rh < 1.)
        result += 15*pow(1-rh,5);

      double sigma = quintic_spline_sigma[gdimension-1]
                   / pow(h,gdimension);
      result *= sigma;
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
    const point_t & vecP,
    const double h)
  {
    const double r = vector_norm(vecP);
    double rh = 3.*r/h;

    point_t result = 0.0;
    if (rh < 3.) {
      double sigma = 3.*quintic_spline_sigma[gdimension-1]
                   / pow(h,gdimension+1);
      double dWdr = -5.*pow(3-rh,4);
      if(rh < 2.)
        dWdr += 30.*pow(2-rh,4);

      if(rh < 1.)
        dWdr += -75.*pow(1-rh,4);

      result = vecP*sigma*dWdr/r;
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
    const double r, 
    const double h)
  {
    double rh = r/h;
    //double rh = 2.*r/h;
    double result = 0.;

    if(rh < 1.0) {
    //if(rh < 2.0) {
      double rh2 = (1 - rh)*(1 - rh);
      //double rh2 = (1 - rh/2.)*(1 - rh/2.);
      double sigma = wendland_quintic_sigma[gdimension-1]
                   / pow(h,gdimension);

      // Different cases for 1D and 2D/3D 
      if(gdimension == 1)
        result = sigma*rh2*(1 - rh)*(3*rh + 1);
      else
        result = sigma*rh2*rh2*(4*rh + 1);
        //result = sigma*rh2*rh2*(2*rh + 1);
    }
    return result; 
  }

  /**
   * @brief      Gradient of Wendland quintic kernel
   * \TODO add the ref
   *
   * @param[in]  vecP  The vector pab = pa - pb 
   * @param[in]  h     The smoothing length 
   *
   * @return     Contribution from the particle 
   */
  point_t 
  gradient_wendland_quintic(
    const point_t & vecP,
    const double h)
  {
    double r = vector_norm(vecP);
    double rh = r/h;
    //double rh = 2.*r/h;
    point_t result = 0.0;

    if(rh < 1.0) {
    //if(rh < 2.0) {
      double rh2 = (1 - rh)*(1 - rh);
      //double rh2 = (1 - rh/2.)*(1 - rh/2.);
      double sigma = 2.*wendland_quintic_sigma[gdimension-1]
                   / pow(h,gdimension+1);
      double dWdr;
      if(gdimension == 1)
        dWdr = -6.*rh*rh2;
      else 
        dWdr = -10.*rh*rh2*(1 - rh);
        //dWdr = rh2*(1 - rh/2.)*(-5.*rh); 

      result = vecP*sigma*dWdr/r; 
    }

    return result; 
  }


}; // kernel

#endif // _physics_kernel_h_
