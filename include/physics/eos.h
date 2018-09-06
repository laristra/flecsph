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
 * @brief Create namespace for both analytic and tabulated EOS
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
   * Default values
   * They are usually modified in the main_driver 
   * to be application-specific. 
   * \TODO add a parameter file to be read in the main_driver
   */
  point_t max_boundary = {};
  point_t min_boundary = {};
  double dt = 0.0;
  double damp = 1;
  double totaltime = 0.0;
  double A = 1.0;
  double MAC = 1.;
  int64_t iteration = 0;

  /**
   * @brief      Compute the pressure
   * Ideal gas EOS
   * @param      srch  The source's body holder
   */
  void 
  compute_pressure(
      body_holder* srch)
  { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = (poly_gamma-1.0)*
      source->getDensity()*source->getInternalenergy();
    source->setPressure(pressure);
  } // compute_pressure

#ifdef ADIABATIC
  /**
   * @brief      Compute the pressure based on adiabatic index
   *
   * @param      srch  The source's body holder
   */
  void 
  compute_pressure_adiabatic(
      body_holder* srch)
  { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = source->getAdiabatic()*
      pow(source->getDensity(),poly_gamma);
    source->setPressure(pressure);
  } // compute_pressure
#endif 

  /**
   * @brief      Zero temperature EOS from Chandrasechkar's 
   * 		 This can be used white dwarf system
   *
   * @param      srch  The srch
   */
  void 
  compute_pressure_wd(
      body_holder* srch)
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

#if 0
  /**
   * @brief      Compute the sound speed
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics 
   * 
   * @param      srch  The source's body holder
   */
  void 
  compute_soundspeed(
      body_holder* srch)
  {
    using namespace param;
    body* source = srch->getBody();
    double soundspeed = sqrt(poly_gamma*source->getPressure()
                                       /source->getDensity());
    source->setSoundspeed(soundspeed);
  } // computeSoundspeed

  /**
   * @brief      Compute the density, EOS and spundspeed in the same function 
   * reduce time to gather the neighbors
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void 
  compute_density_pressure_soundspeed(
    body_holder* srch, 
    std::vector<body_holder*>& nbsh)
  {
    compute_density(srch,nbsh);
    compute_pressure(srch);
    compute_soundspeed(srch); 
  }

#ifdef ADIABATIC
  void 
  compute_density_pressure_adiabatic_soundspeed(
    body_holder* srch, 
    std::vector<body_holder*>& nbsh)
  {
    compute_density(srch,nbsh);
    compute_pressure_adiabatic(srch);
    compute_soundspeed(srch); 
  }
#endif

#endif

#if 0 
  // \TODO VERSION USED IN THE BNS, CHECK VALIDITY REGARDING THE OTHER ONE 
  void 
  leapfrog_integration(
      body_holder* srch)
  {
    body* source = srch->getBody();
    

    point_t velocity = source->getVelocityhalf()+
      source->getAcceleration() * dt / 2.;
    point_t velocityHalf = velocity+
      source->getAcceleration() * dt / 2.;
    point_t position = source->getPosition()+velocityHalf*dt;
    // integrate dadt 
    double adiabatic_factor = source->getAdiabatic() + source->getDadt()* dt;

    source->setVelocity(velocity);
    source->setVelocityhalf(velocityHalf);
    source->setPosition(position);
    source->setAdiabatic(adiabatic_factor);
    
    mpi_assert(!std::isnan(position[0])); 
  }
#endif

}; // eos

#endif // _eos_h_
