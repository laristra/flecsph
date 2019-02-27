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
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5;

    double x_wd = pow((source.getDensity())/B_wd,1.0/3.0);
    double pressure = A_wd*(x_wd*(2.0*x_wd*x_wd-3.0)*
 		      sqrt(x_wd*x_wd+1.0)+3.0*asinh(x_wd));
    source.setPressure(pressure);
  } // compute_pressure_wd

// HL : Compute pressure from tabulated EOS. Linking to static lib somewhat
//      need to be fixing

// EOS prep stage for fill eos info from table

#if 0
  void
  EOS_prep(body& source)
  {
    EOS_SC_fill(source.getDensity(), source.getInternalenergy(),
                source.getElectronfraction(),
                1.0//This field should be field for eos cache);
  //May need different source field?
  }
#endif
  void
  compute_pressure_sc(body& source)
  {
    using namespace param;
    //double pressure = EOS_pressure_rho0_u(source->eoscache());
    //source->setPressure(pressure)
  } // compute_pressure_sc

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
   * @brief      Compute sound speed for wd eos
   *
   * @param      srch  The source's body holder
   */
  void compute_soundspeed_wd(body& source) {
    using namespace param;
    double A_wd = 6.00288e22;
    double B_wd = 9.81011e5;
    double x_wd = pow((source.getDensity())/B_wd,1./3.);

    double numer = 8.*source.getDensity()*x_wd - 3.*B_wd;
    double deno = 3*B_wd*B_wd*x_wd*x_wd*sqrt(x_wd*x_wd+1);

    double soundspeed = A_wd*(numer/deno +
                              x_wd/(3.*source.getDensity()
                                    *sqrt(1-x_wd*x_wd)));
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
  else if(boost::iequals(eos_type, "stellar collapse")) {
    init = init_ideal;  // TODO
    compute_pressure = compute_pressure_sc;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else {
    std::cerr << "Bad eos_type parameter" << std::endl;
  }
}

} // namespace eos

#endif // _eos_h_
