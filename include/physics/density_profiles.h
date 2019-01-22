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
 * @file density_profiles.h
 * @author Oleg Korobkin
 * @date January 2019
 * @brief interface to select various density profiles
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

  /**
   * @brief  parabolic density
   * @param  r     - spherical radius
   */
  double rho_parabolic_density(const double r) {
    return 15./(8.*M_PI)*(1. - r*r);
  }

  double mass_parabolic_density(const double r) {
    return 7.5*CU(r)*(1./3. - r*r/5.);
    return CU(r);
  }

  /**
   * @brief      External force selector
   * @param      efstr    ext. force string
   */
  void select() {
    using namespace param;
    if (boost::iequals(density_profile,"constant")) {
      spherical_density_profile = rho_constant_density;
      spherical_mass_profile = mass_constant_density;
    }
    else if (boost::iequals(density_profile,"parabolic")) {
      spherical_density_profile = rho_parabolic_density;
      spherical_mass_profile = mass_parabolic_density;
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
