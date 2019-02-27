/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@ @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@///// /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@      /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@ /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////  /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@      /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@      /@@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file physics.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic physics implementation
 */

#ifndef _BNS_physics_h_
#define _BNS_physics_h_

#include <vector>

#include "params.h"
#include "kernels.h"
#include "tree.h"

namespace physics{

  // Specific for BNS
  double angular_moment = 1.36049047255;
  double angular_speed = 1.0;
  double QZZ = 0.;
  double Omega2 = 0.;
  double t_relax = 1.12;

  void 
  compute_QZZ(
    std::vector<body_holder*>& bodies)
  {
    QZZ = 0.;
    Omega2 = 0.;
    for(auto bh: bodies){
      body * src = bh->getBody();
      QZZ += src->getMass() * 
        (src->getPosition()[0]*src->getPosition()[0]+
          src->getPosition()[1]*src->getPosition()[1]);
    }
  }

  void compute_rotation(
    std::vector<body_holder*>& bodies)
  {
    Omega2 = pow(angular_moment/QZZ,2);

    for(auto bh: bodies){
      body * src = bh->getBody();
      point_t frot = src->getPosition();
      frot[0] *= Omega2;
      frot[1] *= Omega2;
      if(gdimension == 3){
        frot[2] = 0.;
      }
      src->setAcceleration(src->getAcceleration() + frot);
      // Correct the velocity 
      src->setAcceleration(src->getAcceleration() - 
        src->getVelocity()/t_relax);
    }
  }
  
  void 
  apply_rotation(
    body_holder* srch)
  {
    body* source = srch->getBody();
    double theta = angular_speed * dt;
    mpi_assert(theta > 0);
    point_t tmp = source->getPosition();
    point_t new_position = tmp;
    new_position[0] = tmp[0]*cos(theta)-tmp[1]*sin(theta);
    new_position[1] = tmp[0]*sin(theta)+tmp[1]*cos(theta);
    source->setPosition(new_position);
  }

  // Internal energy from Adiabatic factor 
  // u = A / (g-1) * rho ^ {g - 1}
  void 
  compute_internal_energy_from_adiabatic(
    body_holder* srch)
  {
    using namespace param;
    body* source = srch->getBody();
    double u = source->getAdiabatic()/(poly_gamma-1)*
      pow(source->getDensity(),poly_gamma-1);
    source->setInternalenergy(u);
  }

  // Adiabatic ratio from internal energy 
  void 
  compute_adiabatic_from_internal_energy(
    body_holder* srch)
  {
    using namespace param;
    body* source = srch->getBody();
    double A = (poly_gamma-1)*source->getInternalenergy()/
      pow(source->getDensity(),poly_gamma-1);
    source->setAdiabatic(A);
  }

}; // physics

#endif // _BNS_physics_h_
