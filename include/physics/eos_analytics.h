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
 * @file eos_analytics.h
 * @author Julien Loiseau, Hyun Lim
 * @date June 2017 , Jul 2017 (for more analytic EOS)
 * @brief Implementation of analytics EOS
 */

#ifndef _physics_eos_analytics_h_
#define _physics_eos_analytics_h_

#include <vector>

#include "eos.h"

class eos_analytics:
 public eos{
  
public:
  eos_analytics(){};
  ~eos_analytics(){};
  
  eos_analytics(double gamma): eos(gamma){};

  //static
  double compute_pressure(
    body_holder* srch
  ){
    body* source = srch->getBody(); 
    double pressure = (gamma_-1.0)*
      (source->getDensity())*(source->getInternalenergy());
    source->setPressure(pressure); 
  };
  //Zero temperature EOS for double white dwarf
  double compute_pressure_wd(//HL : Need to merge as one function
    body_holder* srch
  ){
    body* source = srch->getBody();
    double A_dwd = 6.00288e22; 
    double B_dwd = 9.81011e5;

    double x_dwd = pow((source->getDensity())/B_dwd,1.0/3.0); 
    double pressure = A_dwd*(x_dwd*(2.0*x_dwd*x_dwd-3.0)*
                      pow(x_dwd*x_dwd+1.0,1.0/2.0)+3.0*asinh(x_dwd));
    source->setPressure(pressure); 
  };

private: 

};

#endif // _physics_eos_analytics_h_
