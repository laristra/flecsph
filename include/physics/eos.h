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


namespace eos{
  using namespace param;

  /**
   * @brief      Compute the pressure for ideal gas EOS
   * @param      srch  The source's body holder
   */
  void compute_pressure_ideal(body_holder* srch) { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = (poly_gamma-1.0)*
      source->getDensity()*source->getInternalenergy();
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
   * 		 This can be used white dwarf system
   *
   * @param      srch  The srch
   */
  void compute_pressure_wd( body_holder* srch)
  { 
    body* source = srch->getBody();
    double A_dwd = 6.00288e22;
    double B_dwd = 9.81011e5;

    double x_dwd = pow((source->getDensity())/B_dwd,1.0/3.0);
    double pressure = A_dwd*(x_dwd*(2.0*x_dwd*x_dwd-3.0)*
 		      pow(x_dwd*x_dwd+1.0,1.0/2.0)+3.0*asinh(x_dwd));
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

  // TODO: add soundspeed for WDs and tabulates SC eos

  
  // eos function types and pointers
  typedef void (*compute_pressure_t)(body_holder*);
  typedef void (*compute_soundspeed_t)(body_holder*);
  compute_pressure_t compute_pressure = compute_pressure_ideal;
  compute_soundspeed_t compute_soundspeed = compute_soundspeed_ideal;

}; // namespace eos

#endif // _eos_h_
