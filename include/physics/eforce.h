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
  point_t acceleration(body_holder* srch) { 
    // body* source = srch->getBody();
    point_t a = 0.0;
    return a;
  }


  /**
   * @brief      Trivial case: zero potential
   * @param      srch  The source's body holder
   */
  double potential(body_holder* srch)
  { 
    // body* source = srch->getBody();
    return 0.0;
  }


} // namespace external_force

#endif // _eforce_h_
