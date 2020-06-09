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
 * @brief Namespace for both analytic EOS
 *        for pressure computation and sound speed
 */

#pragma once

#include <vector>

#include "params.h"
#include "tree.h"
#include "utils.h"
#include <boost/algorithm/string.hpp>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

namespace eos {
using namespace param;

/**
 * @brief      Lazy function: does nothing
 * @param      source    source's body holder
 */
void
dummy(body & source) {}

/////////////////////////////////////////////////////////////////////////////
// EOS INITIALIZATION

/**
 * @brief      Equation-of-state intializer:
 *             computes missing quantities etc.
 * @param      srch  The source's body holder
 */
void
init_polytropic(body & source) {
  using namespace param;
  double K = source.getPressure() / pow(source.getDensity(), poly_gamma);
  source.setAdiabatic(K);
  return;
}

/////////////////////////////////////////////////////////////////////////////
// PRESSURE

/**
 * @brief      Compute the pressure for ideal gas EOS:
 *             P(\rho, u) = (\Gamma - 1)*\rho*u
 *
 * @param      srch  The source's body holder
 */
void
compute_pressure_ideal(body & source) {
  using namespace param;
  double pressure =
    (poly_gamma - 1.0) * source.getDensity() * source.getInternalenergy();
  source.setPressure(pressure);
}

/**
 * @brief      Compute pressure from density using polytrope
 *             P(\rho) = A*\rho^\Gamma
 *
 * @param      srch  The source's body holder
 */
void
compute_pressure_poly(body & source) {
  using namespace param;
  double pressure =
    source.getAdiabatic() * pow(source.getDensity(), poly_gamma);
  source.setPressure(pressure);
}

/**
 * @brief      Zero-T electron-denegerate EOS (after Chandrasechkar)
 * 		     Used for a cold white dwarf
 *
 * @param      source  The particle
 */
void
compute_pressure_wd(body & source) {
  double Ye = 0.5;
  double A_wd = 6.00288e22;
  double B_wd = 9.81011e5 / Ye;

  double x_wd = pow((source.getDensity()) / B_wd, 1.0 / 3.0);
  double pressure =
    A_wd *
    (x_wd * (2.0 * SQ(x_wd) - 3.0) * sqrt(SQ(x_wd) + 1.0) + 3.0 * asinh(x_wd));
  source.setPressure(pressure);
} // compute_pressure_wd

/**
 * @brief      Compute the pressure for piecewise-polytrope EOS
 * @param      srch  The source's body holder
 */
void
compute_pressure_ppt(body & source) {
  using namespace param;

  // TODO : transient density might be parametrized
  //       or determining automatically.
  //       But here I put certian value that is
  //       reasonalbe transition between
  //       relativisitc and non-relativistic regimes
  const double transition_density = 5e+14;

  if(source.getDensity() <= transition_density) {
    double pressure =
      source.getAdiabatic() * pow(source.getDensity(), poly_gamma);
    source.setPressure(pressure);
  }
  else {
    double pressure =
      (source.getAdiabatic() * pow(transition_density, poly_gamma) /
        pow(transition_density, poly_gamma2)) *
      pow(source.getDensity(), poly_gamma2);
    source.setPressure(pressure);
  }
} // compute_pressure_ppt

/////////////////////////////////////////////////////////////////////////////
// SOUND SPEED

/**
 * @brief      Compute sound speed for ideal fluid or polytropic eos
 * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
 *
 * @param      srch  The source's body holder
 */
void
compute_soundspeed_ideal(body & source) {
  using namespace param;
  double soundspeed =
    sqrt(poly_gamma * source.getPressure() / source.getDensity());
  source.setSoundspeed(soundspeed);
}

/**
 * @brief      Compute sound speed for polytropic eos
 * From "Relativistic Hydrodynamics", Rezzolla & Zanotti
 * TODO: this is never used
 *
 * @param      srch  The source's body holder
 */
void
compute_soundspeed_poly(body & source) {
  using namespace param;
  const double cc = 2.99792458e10; // [cm/s] TODO
  const double P = source.getPressure(), rho = source.getDensity(),
               u = source.getInternalenergy();
  double enth = SQ(cc) + u + P / rho;
  source.setSoundspeed(sqrt(poly_gamma * P / (rho * enth)) * cc);
} // compute_soundspeed_poly

/**
 * @brief      Compute sound speed for piecewise polytropic eos
 *
 * @param      srch  The source's body holder
 */
void
compute_soundspeed_ppt(body & source) {
  using namespace param;
  double density_transition = 500000000000000;
  if(source.getDensity() <= density_transition) {
    double soundspeed =
      sqrt(poly_gamma * source.getPressure() / source.getDensity());
    source.setSoundspeed(soundspeed);
  }
  else {
    double soundspeed =
      sqrt(poly_gamma2 * source.getPressure() / source.getDensity());
    source.setSoundspeed(soundspeed);
  }
} // compute_soundspeed_ppt

/**
 * @brief      Compute sound speed for wd eos
 *
 * @param      srch  The source's body holder
 */
void
compute_soundspeed_wd(body & source) {
  using namespace param;
  double Ye = 0.5;
  double A_wd = 6.00288e22;
  double B_wd = 9.81011e5 / Ye;
  double x_wd = pow((source.getDensity()) / B_wd, 1. / 3.);
  double cc = 29979245800.0;

  double sterm = sqrt(1.0 + SQ(x_wd));
  double numer = 3.0 / sterm + sterm * (6.0 * SQ(x_wd) - 3.0) +
                 SQ(x_wd) * (2.0 * SQ(x_wd) - 3.0) / sterm;
  double denom = -1.0 / sterm + sterm * (1.0 + 6.0 * SQ(x_wd)) +
                 SQ(x_wd) * (1.0 + 2.0 * SQ(x_wd)) / sterm;

  double soundspeed = sqrt(numer / (3.0 * denom)) * cc;

  if(not(numer / denom > 0)) {
    std::cout << "speed of sounds is not a real number: "
              << "numer/denom = " << numer / denom << std::endl;
    std::cout << "Failed particle id: " << source.id() << std::endl;
    std::cerr << "particle position: " << source.coordinates() << std::endl;
    std::cerr << "particle velocity: " << source.getVelocity() << std::endl;
    std::cerr << "particle acceleration: " << source.getAcceleration()
              << std::endl;
    std::cerr << "smoothing length:  " << source.radius() << std::endl;
    assert(false);
  }

  // double numer = 8.*source.getDensity()*x_wd - 3.*B_wd;
  // double deno = 3*B_wd*B_wd*x_wd*x_wd*sqrt(x_wd*x_wd+1);

  // double soundspeed = A_wd*(numer/deno +
  //                          x_wd/(3.*source.getDensity()
  //                                *sqrt(1-x_wd*x_wd)));
  source.setSoundspeed(soundspeed);
} // compute_soundspeed_wd

/////////////////////////////////////////////////////////////////////////////
// TEMPERATURE

/**
 * @brief      Compute temperature via ideal gas in C/O WD
 *             TODO: parameterize abar, zbar; double-check formula [???]
 *
 * @param      srch  The source's body holder
 */
void
compute_temperature_idealgas(body & source) {
  const double kB = 1.3806505e-16, // [erg/K]
    abar = 12.0, // [mol/g] molar mass of Carbon-12
    zbar = 6.0, // proton number for C
    amu = 1.66053906660e-24, // [g] a.m.u.
    me = 9.10938356e-28; // [g] electron mass
  const double P = source.getPressure(), rho = source.getDensity(),
               Ye = source.getElectronfraction();
  double mu = abar * (amu + Ye * me) / (zbar + 1.0); // ???
  double T = mu * P / (rho * kB);
  source.setTemperature(T);
} // compute_temperature_ideal

/////////////////////////////////////////////////////////////////////////////
// EOS SELECTORS

// eos function types and pointers
typedef void (*compute_quantity_t)(body &);
compute_quantity_t compute_pressure = compute_pressure_ideal;
compute_quantity_t compute_soundspeed = compute_soundspeed_ideal;
compute_quantity_t compute_temperature = dummy;

typedef void (*eos_init_t)(body &);
eos_init_t init = dummy;

/**
 * @brief  Installs the 'compute_pressure' and 'compute_soundspeed'
 *         function pointers, depending on the value of eos_type
 */
void
select(const std::string & eos_type) {
  if(boost::iequals(eos_type, "ideal fluid")) {
    init = dummy;
    compute_pressure = compute_pressure_ideal;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else if(boost::iequals(eos_type, "polytropic")) {
    init = init_polytropic;
    compute_pressure = compute_pressure_poly;
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else if(boost::iequals(eos_type, "white dwarf")) {
    init = dummy; // TODO
    compute_pressure = compute_pressure_wd;
    compute_soundspeed = compute_soundspeed_wd;
    compute_temperature = compute_temperature_idealgas;
  }
  else if(boost::iequals(eos_type, "piecewise polytropic")) {
    init = dummy; // TODO
    compute_pressure = compute_pressure_ppt;
    compute_soundspeed = compute_soundspeed_ppt;
  }
  else if (boost::iequals(eos_type, "no eos")) {
    // This eos does nothing
    init = dummy;
    compute_pressure = dummy;
    compute_soundspeed = dummy;
  }
  else if(boost::iequals(eos_type, "pure gravitational collapse")) {
    // This eos only calcultes the soundspeed
    // and the pressure is kept at the initial value
    // It's for pure collapse simulations so that pressure does not
    // counteract the collapse
    init = dummy;
    compute_pressure = dummy; // do nothing to the pressure
    compute_soundspeed = compute_soundspeed_ideal;
  }
  else {
    std::cerr << "Bad eos_type parameter" << std::endl;
  }
}

} // namespace eos
