/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file density_profiles.h
 * @author Oleg Korobkin
 * @date January 2019
 * @brief Interface to select various density profiles
 * 
 * Notes: 
 *  - all spherical density profiles have support with radius R = 1;
 *  - all profiles are normalized to total mass M = 1.
 */

#ifndef DENSITY_PROFILES_H
#define DENSITY_PROFILES_H

#include <stdlib.h>
#include "user.h"
#include "tree.h"
#include <math.h>
#include <boost/algorithm/string.hpp>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))
namespace density_profiles {

  // spherical density profile function
  typedef double  (*radial_function_t)(const double);
  static radial_function_t spherical_density_profile = NULL;
  static radial_function_t spherical_mass_profile = NULL;
  static radial_function_t spherical_drho_dr = NULL;

  // an important constant for the mesa density
  static double mesa_rho0;
  const  double mesa_q = 0.25; // ratio of the slope width to the radius

  /**
   * @brief  constant uniform density in a domain of radius R = 1
   * @param  r     - spherical radius
   */
  double rho_constant_density(const double r) {
    return 0.75/M_PI;
  }

  double mass_constant_density(const double r) {
    return CU(r);
  }

  double drhodr_constant_density(const double r) {
    return 0.;
  }

  /**
   * @brief  parabolic density
   * @param  r     - spherical radius
   */
  double rho_parabolic_density(const double r) {
    return 15./(8.*M_PI)*(1. - r*r);
  }

  double mass_parabolic_density(const double r) {
    return 7.5*CU(r)*(1./3. - r*r/5.);
  }

  double drhodr_parabolic_density(const double r) {
    return -15./(4.*M_PI)*r;
  }

  /**
   * @brief  spherical "mesa" density: flat top and steep slopes
   *
   *           / rho0                      if r < r0;
   *           |
   * rho(r) = <  rho0 (1 - (r-r0)^2/dr^2)  if r < r0 + dr;
   *           |
   *           \ 0                         if r > r0 + dr.
   *
   * @param  r     - spherical radius
   */
  double mesa_mass_helper(const double r) { 
    const double dr = mesa_q, r0 = 1. - mesa_q;
    return CU(r0)/3.
       + (CU(r)-CU(r0))/3.*(1. - SQ(r0)/SQ(dr))
       + r0*(SQ(r)-SQ(r0))*(SQ(r)+SQ(r0))/(2.*SQ(dr))
       - (SQ(r)*CU(r)-SQ(r0)*CU(r0))/(5*SQ(dr)); 
  }

  double rho_mesa_density(const double r) {
    double rho = 0.0;
    const double r0 = 1. - mesa_q;
    if (r < 1.-mesa_q) 
      rho = mesa_rho0;
    else if (r < 1.)
      rho = mesa_rho0*(1. - SQ(r-r0)/SQ(mesa_q));
    return rho;
  }

  double mass_mesa_density(const double r) {
    double m = 0.0;
    if (r < 1.-mesa_q) 
      m = 4.*M_PI/3. * mesa_rho0 * CU(r);
    else if (r - 1. < 1e-12)
      m = 4.*M_PI*mesa_rho0*mesa_mass_helper(r);
    return m;
  }

  double drhodr_mesa_density(const double r) {
    double drhodr = 0.0;
    const double r0 = 1. - mesa_q;
    if (r > 1.-mesa_q and r < 1.) 
      drhodr = -2.*mesa_rho0*(r-r0)/SQ(mesa_q);
    return drhodr;
  }

  /**
   * @brief      Density profile selector
   */
  void select() {
    using namespace param;
    if (boost::iequals(density_profile,"constant")) {
      spherical_density_profile = rho_constant_density;
      spherical_mass_profile = mass_constant_density;
      spherical_drho_dr = drhodr_constant_density;
    }
    else if (boost::iequals(density_profile,"parabolic")) {
      spherical_density_profile = rho_parabolic_density;
      spherical_mass_profile = mass_parabolic_density;
      spherical_drho_dr = drhodr_parabolic_density;
    }
    else if (boost::iequals(density_profile,"mesa")) {
      spherical_density_profile = rho_mesa_density;
      spherical_mass_profile = mass_mesa_density;
      spherical_drho_dr = drhodr_mesa_density;
      mesa_rho0 = 1./(4.*M_PI*mesa_mass_helper(1.));
    }
    else {
      clog(error) << "ERROR: wrong parameter in density_profiles";
      exit(2);
    }

  } // select()

} // namespace density_profiles

#undef SQ
#undef CU
#endif // DENSITY_PROFILES_H
