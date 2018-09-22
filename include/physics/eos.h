/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
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
 *        for pressure computation
 */

#ifndef _eos_h_
#define _eos_h_

#include <vector>

#include "params.h"
#include "utils.h"
#include "tree.h"
#include <boost/algorithm/string.hpp>


namespace eos {
  using namespace param;

  /**
   * @brief      Compute the pressure for ideal gas EOS
   * @param      srch  The source's body holder
   */
  void compute_pressure_ideal(body_holder* srch) { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = (poly_gamma-1.0)*source->getDensity() 
                                      *source->getInternalenergy();
    source->setPressure(pressure);
  }


  /**
   * @brief      Compute the pressure based on adiabatic index
   * @param      srch  The source's body holder
   */
  void compute_pressure_adiabatic( body_holder* srch)
  { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = source->getAdiabatic()*
      pow(source->getDensity(),poly_gamma);
    source->setPressure(pressure);
  }

  /**
   * @brief      Zero temperature EOS from Chandrasechkar's 
   * 		     This can be used for a cold white dwarf
   *
   * @param      srch  The srch
   */
  void compute_pressure_wd( body_holder* srch)
  { 
    body* source = srch->getBody();
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5;

    double x_wd = pow((source->getDensity())/B_wd,1.0/3.0);
    double pressure = A_wd*(x_wd*(2.0*x_wd*x_wd-3.0)*
 		      sqrt(x_wd*x_wd+1.0)+3.0*asinh(x_wd));
    source->setPressure(pressure);
  } // compute_pressure_wd

// HL : Compute pressure from tabulated EOS. Working now..

#if 0

#ifdef _EOS_TAB_SC

  void
  EOS_prep(body_holder* srch)
  {
    body* source = srch->getBody();
    EOS_SC_Fill();
  }

  void
  compute_pressure_tabEOS_SC(
      body_holder* srch)
  { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = EOS_pressure();
  } // compute_pressure_tabEOS_SC

#endif

#endif

  /**
   * @brief      Compute sound speed for ideal fluid or polytropic eos
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics 
   * 
   * @param      srch  The source's body holder
   */
  void compute_soundspeed_ideal(body_holder* srch) {
    using namespace param;
    body* source = srch->getBody();
    double soundspeed = sqrt(poly_gamma*source->getPressure()
                                       /source->getDensity());
    source->setSoundspeed(soundspeed);
  }

  /**
   * @brief      Compute sound speed for wd eos
   * 
   * @param      srch  The source's body holder
   */
  void compute_soundspeed_wd(body_holder* srch) {
    using namespace param;
    body* source = srch->getBody();
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5;
    double x_wd = pow((source->getDensity())/B_wd,1./3.);

    double numer = 8.*source->getDensity()*x_wd - 3.*B_wd;
    double deno = 3*B_wd*B_wd*x_wd*x_wd*sqrt(x_wd*x_wd+1);

    double soundspeed = A_wd*(numer/deno + 
                              x_wd/(3.*source->getDensity()
                                    *sqrt(1-x_wd*x_wd)));
    source->setSoundspeed(soundspeed);
  }

  // TODO: add tabulated eos from SC
  
  // eos function types and pointers
  typedef void (*compute_pressure_t)(body_holder*);
  typedef void (*compute_soundspeed_t)(body_holder*);
  compute_pressure_t compute_pressure = compute_pressure_ideal;
  compute_soundspeed_t compute_soundspeed = compute_soundspeed_ideal;

/**
 * @brief  Installs the 'compute_pressure' and 'compute_soundspeed' 
 *         function pointers, depending on the value of eos_type
 */
void select(const std::string& eos_type) {
  if(boost::iequals(eos_type, "ideal fluid")) {
    compute_pressure = compute_pressure_ideal;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else if(boost::iequals(eos_type, "polytropic")) {
    compute_pressure = compute_pressure_adiabatic;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else if(boost::iequals(eos_type, "white dwarf")) {
    compute_pressure = compute_pressure_wd;
    compute_soundspeed = compute_soundspeed_wd;
  }
  else {
    std::cerr << "Bad eos_type parameter" << std::endl;
  }
}

} // namespace eos

#endif // _eos_h_
