/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
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
 * @file eos.h
 * @brief Namespace for both analytic and tabulated EOS
 *        for pressure computation and sound speed
 */

#ifndef _eos_h_
#define _eos_h_

#include <vector>

#include "params.h"
#include "utils.h"
#include "tree.h"
#include <boost/algorithm/string.hpp>

// HL : This is for adding tabulated EOS redaer
//#include "eos_stellar_collapse.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

namespace eos {
  using namespace param;

  /**
   * @brief      Equation-of-state intializer:
   *             computes missing quantities etc.
   * @param      srch  The source's body holder
   */
  void init_ideal(body& source) { return; } // do nothing
  void init_polytropic(body& source) {
    using namespace param;
    double K = source.getPressure()
             / pow(source.getDensity(),poly_gamma);
    source.setAdiabatic(K);
    return;
  }


  /**
   * @brief      Compute the pressure for ideal gas EOS
   * @param      srch  The source's body holder
   */
  void compute_pressure_ideal(body& source) {
    using namespace param;
    double pressure = (poly_gamma-1.0)*source.getDensity()
                                      *source.getInternalenergy();
    source.setPressure(pressure);
  }

  /**
   * @brief      Output pressure same as input pressure
   * @param      srch  The source's body holder
   */
  void compute_pressure_no_eos(body& source) {
    using namespace param;
    double pressure = source.getPressure();
    source.setPressure(pressure);
  }

  /**
   * @brief      Compute the pressure based on adiabatic index
   * @param      srch  The source's body holder
   */
  void compute_pressure_adiabatic( body& source)
  {
    using namespace param;
    double pressure = source.getAdiabatic()*
      pow(source.getDensity(),poly_gamma);
    source.setPressure(pressure);
  }

  /**
   * @brief      Zero temperature EOS from Chandrasechkar's
   * 		     This can be used for a cold white dwarf
   *
   * @param      srch  The srch
   */
  void compute_pressure_wd(body& source)
  {
    double Ye   = 0.5;
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5/Ye;

    double x_wd = pow((source.getDensity())/B_wd,1.0/3.0);
    double pressure = A_wd*(x_wd*(2.0*SQ(x_wd) - 3.0)*
 		      sqrt(SQ(x_wd) + 1.0) + 3.0*asinh(x_wd));
    source.setPressure(pressure);
  } // compute_pressure_wd

  #if 1
  // HL : since we are merging tab EOS, I am adding ppt anyway
  /**
   * @brief      Compute the pressure based on piecewise polytrope
   * @param      srch  The source's body holder
   */
  void compute_pressure_ppt( body& source)
  {
    using namespace param;

    //TODO : transient density might be parametrized
    //       or determining automatically.
    //       But here I put certian value that is
    //       reasonalbe transition between
    //       relativisitc and non-relativistic regimes
    double density_transition = 500000000000000;

    if(source.getDensity() <= density_transition) {
       double pressure = source.getAdiabatic()*
       pow(source.getDensity(),poly_gamma);
       source.setPressure(pressure);
    }
    else {
       double pressure = (source.getAdiabatic()*
       pow(density_transition,poly_gamma)/
       pow(density_transition,poly_gamma2))*
       pow(source.getDensity(),poly_gamma2);
       source.setPressure(pressure);
    }
  } //compute_pressure_ppt
  #endif

#if 0
/************************************************************************/
//May.30.2019
// Start SC EOS reader merging
// This is a pusedo-code. This just shows guideline how we can
// use SC reader to get P and Cs

// EOS prep stage for fill eos info from table
  void
  EOS_prep(body& source)
  {
    EOS_pressure_rho0_u(source);
    EOS_sound_speed_rho0_u(source);
  }

// Getting pressure
  void
  compute_pressure_sc(body& source)
  {
    EOS_pressure_rho0_u(source);
  } // compute_pressure_sc

// Getting soundspeed
  void
  compute_soundspeed_sc(body& source)
  {
    EOS_sound_speed_rho0_u(source);
  }
/***************************************************************************/
#endif

 /**
   * @brief      Compute sound speed for ideal fluid or polytropic eos
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
   *
   * @param      srch  The source's body holder
   */
  void compute_soundspeed_ideal(body& source) {
    using namespace param;
    double soundspeed = sqrt(poly_gamma*source.getPressure()
                                       /source.getDensity());
    source.setSoundspeed(soundspeed);
  }

 /**
   * @brief      Compute sound speed for piecewise polytropic eos
   *
   * @param      srch  The source's body holder
   */
  void compute_soundspeed_ppt(body& source) {
    using namespace param;
    double density_transition = 500000000000000;
    if(source.getDensity() <= density_transition) {
       double soundspeed = sqrt(poly_gamma*source.getPressure()
                                         /source.getDensity());
       source.setSoundspeed(soundspeed);
    }
    else {
       double soundspeed = sqrt(poly_gamma2*source.getPressure()
                                          /source.getDensity());
       source.setSoundspeed(soundspeed);
    }
  }

  /**
   * @brief      Compute sound speed for polytropic eos
   * From "Relativistic Hydrodynamics", Rezzolla & Zanotti
   *
   * @param      srch  The source's body holder
   */
  void compute_soundspeed_adiabatic(body& source) {
    using namespace param;
    double cc = 29979245800.0;
    double enth = cc*cc + source.getInternalenergy() + source.getPressure()/source.getDensity();
    double soundspeed = sqrt(poly_gamma*source.getPressure()/(source.getDensity()*enth))*cc;
    source.setSoundspeed(soundspeed);
  }



  /**
   * @brief      Compute sound speed for wd eos
   *
   * @param      srch  The source's body holder
   */
  void compute_soundspeed_wd(body& source) {
    using namespace param;
    double Ye   = 0.5;
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5/Ye;
    double x_wd = pow((source.getDensity())/B_wd,1./3.);
    double cc   = 29979245800.0;

    double sterm = sqrt(1.0 + SQ(x_wd));
    double numer = 3.0/sterm + sterm*(6.0*SQ(x_wd) - 3.0)
                 + SQ(x_wd) * (2.0*SQ(x_wd) - 3.0)/sterm;
    double denom = -1.0/sterm + sterm*(1.0 + 6.0*SQ(x_wd))
                 + SQ(x_wd)*(1.0 + 2.0*SQ(x_wd))/sterm;

    double soundspeed = sqrt(numer/(3.0*denom))*cc;

    if (not (numer/denom>0)) {
      std::cout << "speed of sounds is not a real number: "
                << "numer/denom = " << numer/denom << std::endl;
      std::cout << "Failed particle id: " << source.id() << std::endl;
      std::cerr << "particle position: " << source.coordinates() << std::endl;
      std::cerr << "particle velocity: " << source.getVelocity() << std::endl;
      std::cerr << "particle acceleration: " << source.getAcceleration() << std::endl;
      std::cerr << "smoothing length:  " << source.radius()
                                         << std::endl;
      assert (false);
    }

    //double numer = 8.*source.getDensity()*x_wd - 3.*B_wd;
    //double deno = 3*B_wd*B_wd*x_wd*x_wd*sqrt(x_wd*x_wd+1);

    //double soundspeed = A_wd*(numer/deno +
    //                          x_wd/(3.*source.getDensity()
    //                                *sqrt(1-x_wd*x_wd)));
    source.setSoundspeed(soundspeed);
  }

  // eos function types and pointers
  typedef void (*compute_quantity_t)(body&);
  compute_quantity_t compute_pressure = compute_pressure_ideal;
  compute_quantity_t compute_soundspeed = compute_soundspeed_ideal;

  typedef void (*eos_init_t)(body&);
  eos_init_t init = init_ideal;

/**
 * @brief  Installs the 'compute_pressure' and 'compute_soundspeed'
 *         function pointers, depending on the value of eos_type
 */
void select(const std::string& eos_type) {
  if(boost::iequals(eos_type, "ideal fluid")) {
    init = init_ideal;
    compute_pressure = compute_pressure_ideal;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else if(boost::iequals(eos_type, "no eos")) {
    init = init_polytropic;
    compute_pressure = compute_pressure_no_eos;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else if(boost::iequals(eos_type, "polytropic")) {
    init = init_polytropic;
    compute_pressure = compute_pressure_adiabatic;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else if(boost::iequals(eos_type, "white dwarf")) {
    init = init_ideal;  // TODO
    compute_pressure = compute_pressure_wd;
    compute_soundspeed = compute_soundspeed_wd;
  }
  else if(boost::iequals(eos_type, "piecewise polytropic")) {
    init = init_ideal;  // TODO
    compute_pressure = compute_pressure_ppt;
    compute_soundspeed = compute_soundspeed_ppt;
  }
  # if 0
  else if(boost::iequals(eos_type, "stellar collapse")) {
    // Reading the table
    init_EOS();
    // Initializing the particles
    init = EOS_prep;  // TODO
    compute_pressure = compute_pressure_sc;
    compute_soundspeed = compute_soundspeed_sc;
  }
  #endif
  else {
    std::cerr << "Bad eos_type parameter" << std::endl;
  }
}

} // namespace eos

#endif // _eos_h_
