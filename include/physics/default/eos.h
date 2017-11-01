/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
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
 * @author Julien Loiseau
 * @date June 2017
 * @brief Abstract EOS implementation
 */

#ifndef _physics_eos_h_
#define _physics_eos_h_

#define DEFAULT_GAMMA 1.4

#include <vector>

#include "tree.h"

class eos{

public:
  eos(double gamma):gamma_(gamma){};
  eos():gamma_(DEFAULT_GAMMA){};

  ~eos(){};

  // Generic eos function to compute the pressure 
  virtual double compute_pressure(
      body_holder*, 
      std::vector<body_holder*>&); 

  virtual double compute_pressure_wd(
      body_holder*, 
      std::vector<body_holder*>&); 

protected:
  double gamma_;
};

#endif // _physics_kernel_h_
