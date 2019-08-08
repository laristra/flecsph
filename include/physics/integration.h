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
#include "utils.h"
#include "user.h"
#include "kernels.h"
#include "tree.h"


namespace integration{
  using namespace param;

  /** 
   * @brief:    Converts cartesian particle coordinates into 
   *            spherical 
   *
   * @param     pos_c   Particle position in cartesian coord.
   */
  point_t
  cartesian_to_spherical (const point_t& pos_c) {
    double r = norm2(pos_c);
    point_t pos_s = 0.0;

    if (gdimension == 1) {
      pos_s[0] = r;
    }
    else if (gdimension == 2) {
      pos_s[0] = r;
      pos_s[1] = atan2(pos_c[1],pos_c[0]);
    }
    else {
      pos_s[0] = r;
      pos_s[1] = atan2(pos_c[1],pos_c[0]);
      pos_s[2] = acos(pos_c[2]/r);
    }
    return pos_s;
  }


  /** 
   * @brief:    Converts spherical particle coordinates into 
   *            cartesian 
   *
   * @param     pos_s   Particle position in spherical coord.
   */
  point_t
  spherical_to_cartesian (const point_t& pos_s) {
    double r = pos_s[0];
    point_t pos_c = 0.0;

    if (gdimension == 1) {
      pos_c[0] = r;
    }
    else if (gdimension == 2) {
      pos_c[0] = r*cos(pos_s[1]);
      pos_c[1] = r*sin(pos_s[1]);
    }
    else {
      pos_c[0] = r*sin(pos_s[2])*cos(pos_s[1]);
      pos_c[1] = r*sin(pos_s[2])*sin(pos_s[1]);
      pos_c[2] = r*cos(pos_s[2]);
    }
    return pos_c;
  }


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

}; //integration 


 #endif // _integration_h_
