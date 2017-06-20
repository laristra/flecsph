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
 * @author Julien Loiseau
 * @date June 2017
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

private: 

};

#endif // _physics_eos_analytics_h_
