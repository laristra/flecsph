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
 * @file integration.h
 * @author Julien Loiseau
 * @date October 2018
 * @brief Integration methods
 */

#ifndef _integration_h_
#define _integration_h_

#include <vector>

#include "params.h"

namespace integration{
  using namespace param;

  /**
   * @brief      Integrate the internal energy variation, update internal energy
   *
   * @param      srch  The source's body holder
   */
  void dadt_integration(
      body& source)
  {
    source.setAdiabatic(
      source.getAdiabatic()+physics::dt*source.getDadt());
  }

  /**
   * @brief      v -> v12
   *
   * @param      srch  The source's body holder
   */
  void
  save_velocityhalf (body& source) {
    source.setVelocityhalf(source.getVelocity());
  }

  /**
   * @brief      Leapfrog: kick velocity
   *             v^{n+1/2} = v^{n} + (dv/dt)^n * dt/2
   *             or
   *             v^{n+1} = v^{n+1/2} + (dv/dt)^n * dt/2
   *
   * @param      srch  The source's body holder
   */
  void
  leapfrog_kick_v (body& source) {
    source.setVelocity(source.getVelocity()
               + 0.5*physics::dt*source.getAcceleration());
  }


  /**
   * @brief      Leapfrog: kick internal energy
   *             u^{n+1/2} = u^{n} + (du/dt)^n * dt/2
   *             or
   *             u^{n+1} = u^{n+1/2} + (du/dt)^n * dt/2
   *
   * @param      srch  The source's body holder
   */
  void
  leapfrog_kick_u (body& source) {
    source.setInternalenergy(source.getInternalenergy()
                     + 0.5*physics::dt*source.getDudt());
  }


  /**
   * @brief      Leapfrog: kick thermokinetic or total energy
   *             e^{n+1/2} = e^{n} + (de/dt)^n * dt/2
   *             or
   *             e^{n+1} = e^{n+1/2} + (de/dt)^n * dt/2
   *
   * @param      srch  The source's body holder
   */
  void
  leapfrog_kick_e (body& source) {
    source.setTotalenergy(source.getTotalenergy()
                     + 0.5*physics::dt*source.getDedt());
  }


  /**
   * @brief      Leapfrog: drift
   *             r^{n+1} = r^{n} + v^{n+1/2} * dt
   *
   * @param      srch  The source's body holder
   */
  void
  leapfrog_drift (body& source) {
    source.set_coordinates(source.coordinates()
                   + physics::dt*source.getVelocity());
  }

}; // integration

#endif // _integration_h_
