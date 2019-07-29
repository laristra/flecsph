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
    point_t rp = source.coordinates() + 0.01*source.radius()*source.getAcceleration();
  
    double r = sqrt(rp[0]*rp[0] + rp[1]*rp[1]);
    double theta = atan2(rp[1],rp[0]); 

    point_t vp = source.getVelocity();
    double vr = cos(theta)*vp[0] + sin(theta)*vp[1];
    double vt = -sin(theta)*vp[0] + cos(theta)*vp[1];

/*
    if (r > 1.0) {
        rp = source.coordinates();
        r = sqrt(rp[0]*rp[0] + rp[1]*rp[1]);
        vp = 0.0;
        source.setVelocity(vp);
    }
*/

    if (r/sphere_radius >= 1.00) {
        if ((rp[0]*vp[0] + rp[1]*vp[1]) <= 0.0) {
            vr = -1.0*vr;
        }
        vp[0] = cos(theta)*vr - r*sin(theta)*vt;
        vp[1] = sin(theta)*vr + r*cos(theta)*vt;
        r = r - 2.0*(r - sphere_radius);
        rp[0] = r*cos(theta);
        rp[1] = r*sin(theta);
        source.setVelocity(vp);
    }

    double rho = (1.0/sqrt(2.0*M_PI*0.1))*(exp(-r*r/(2*0.1)));
    if (r > 1.0) {
      rho = (1.0/sqrt(2.0*M_PI*0.1))*(exp(-1.0*1.0/(2*0.1)));
    }

    double mass = source.getMass();
    double h_a = sph_eta * kernels::kernel_width
          * pow(mass/rho,1./gdimension);
    source.set_radius(h_a);
    source.set_coordinates(rp);
  }

}; // integration

#endif // _integration_h_
