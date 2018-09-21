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
 * @file eforce.h
 * @brief Namespace for the choice of external force and external potential
 */

#ifndef _eforce_h_
#define _eforce_h_

#include <boost/algorithm/string.hpp>

#include "tree.h"
#include "params.h"

namespace external_force {


  /**
   * @brief      Trivial case: zero acceleration
   * @param      srch  The source's body holder
   */
  point_t acceleration_zero(body_holder* srch) { 
    // body* source = srch->getBody();
    point_t a = 0.0;
    return a;
  }

  double potential_zero(body_holder* srch)
  { 
    // body* source = srch->getBody();
    return param::zero_potential_poison_value; // POISON IT
  }

  /**
   * @brief      Square well in y- and z-dimension
   * @param      srch  The source's body holder
   */
  point_t acceleration_squarewell_yz(body_holder* srch) { 
    using namespace param;
    point_t a = 0.0;
    point_t rp =  srch->getBody()->getPosition();
    double box[3];
    box[0] = 0.0;
    box[1] = 0.5*box_width;
    box[2] = 0.5*box_height;
    const double pw_n = extforce_sqwell_power;
    const double pw_a = extforce_sqwell_steepness;
    for (unsigned short i=1; i<gdimension; ++i) {
      a[i]  = (((rp[i] <- box[i]) ? pow(-rp[i]- box[i], pw_n - 1) : 0.0)
              -((rp[i] >  box[i]) ? pow(rp[i] - box[i], pw_n - 1) : 0.0))
            * pw_n*pw_a;
    }
    return a;
  }

  double potential_squarewell_yz(body_holder* srch)
  { 
    // body* source = srch->getBody();
    using namespace param;
    double phi = 0.0;
    point_t rp =  srch->getBody()->getPosition();
    double box[3];
    box[0] = 0.0;
    box[1] = 0.5*box_width;
    box[2] = 0.5*box_height;
    const double pw_n = extforce_sqwell_power;
    const double pw_a = extforce_sqwell_steepness;
    for (unsigned short i=1; i<gdimension; ++i) {
      phi += (((rp[i] <- box[i]) ? pow(-rp[i]- box[i], pw_n) : 0.0)
             +((rp[i] >  box[i]) ? pow(rp[i] - box[i], pw_n) : 0.0))
            *pw_a;
    }
    return phi;
  }


  // acceleration and potential function types and pointers
  typedef double  (*potential_t)(body_holder*);
  typedef point_t (*acceleration_t)(body_holder*);
  potential_t    potential = potential_zero;
  acceleration_t acceleration = acceleration_zero;

  /**
   * @brief      External force selector
   * @param      efstr    ext. force string
   */
  void select(const std::string& efstr) {
    if (boost::iequals(efstr,"zero") or boost::iequals(efstr,"none")) {
      potential = potential_zero;
      acceleration = acceleration_zero;
    }
    else if (boost::iequals(efstr,"square yz-well")) {
      potential = potential_squarewell_yz;
      acceleration = acceleration_squarewell_yz;
    }
    else {
      clog_one(fatal) << "ERROR: bad external_force_type" << std::endl;
    }
  }

  /**
   * @brief      add external potential to the internal energy
   * @param      srch  The source's body holder
   */
  void adjust_internal_energy (body_holder* srch) {
    body* source = srch->getBody();
    source->setInternalenergy(source->getInternalenergy() + potential(srch));
  }

} // namespace external_force

#endif // _eforce_h_
