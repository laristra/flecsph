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

#include <vector>

#include "params.h"
#include "utils.h"
#include "tree.h"


namespace external_force {

  /**
   * @brief      Trivial case: zero acceleration
   * @param      srch  The source's body holder
   */
  point_t acceleration_zero(body_holder* srch) { 
    // body* source = srch->getBody();
    point_t a = {};
    a = 0.0;
    return a;
  }


  /**
   * @brief      Trivial case: zero potential
   * @param      srch  The source's body holder
   */
  double potential_zero(body_holder* srch)
  { 
    // body* source = srch->getBody();
    return -666.0; // POISON IT
  }

  // acceleration and potential function types and pointers
  typedef double  (*potential_t)(body_holder*);
  typedef point_t (*acceleration_t)(body_holder*);
  potential_t    potential = potential_zero;
  acceleration_t acceleration = acceleration_zero;

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
