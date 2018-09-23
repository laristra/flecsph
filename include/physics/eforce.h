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
#define SQ(x) ((x)*(x))

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

  double potential_zero(body_holder* srch) {
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

  double potential_squarewell_yz(body_holder* srch) {
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


  /**
   * @brief      2D airfoil in a hydrodynamic tube
   * @param      srch  The source's body holder
   */
  point_t acceleration_airfoil (body_holder* srch) {
    using namespace param;
    point_t a = 0.0;
    assert (gdimension > 1);
    
    point_t rp =  srch->getBody()->getPosition();
    const double x1 = rp[0] + 1.6;
    const double y1 = rp[1] + 0.3;
    const double alpha = 10*M_PI/180.0;
    const double x = x1*cos(alpha) + y1*sin(alpha);
    const double y =-x1*sin(alpha) + y1*cos(alpha);
    const double pw_n = extforce_sqwell_power;
    const double pw_a = extforce_sqwell_steepness;

    bool inside_bounding_box = std::abs(y)<0.25 && x>-0.05 && x<2.05;
    double upper_surface = .05*x*sqrt(4 - x*x);
    double camber_line   = .1*sin(M_PI*x/2.);
    double phi = SQ(upper_surface) - SQ(y-camber_line) + 0.002;
    if (inside_bounding_box && phi>0.0) {
      double a0, a1;
      a0 = pw_n*pw_a*pow(phi,pw_n-1)
           * ( 2.*(y-camber_line)*(-.1*M_PI/2.*cos(M_PI/2.*x))
             - .0025*(4*x*(2 - x*x)));
      a1 = pw_n*pw_a*pow(phi,pw_n-1) * 2.*(y-camber_line);
      a[0] = a0*cos(alpha) - a1*sin(alpha);
      a[1] = a0*sin(alpha) + a1*cos(alpha);
    }
    return acceleration_squarewell_yz(srch) + a;
  }

  double potential_airfoil (body_holder* srch)
  {
    using namespace param;
    double phi = 0.0;
    assert (gdimension > 1);
    
    point_t rp =  srch->getBody()->getPosition();
    const double x1 = rp[0] + 1.6;
    const double y1 = rp[1] + 0.3;
    const double alpha = 10*M_PI/180.0;
    const double x = x1*cos(alpha) + y1*sin(alpha);
    const double y =-x1*sin(alpha) + y1*cos(alpha);
    const double pw_n = extforce_sqwell_power;
    const double pw_a = extforce_sqwell_steepness;

    bool inside_bounding_box = std::abs(y)<0.25 && x>-0.05 && x<2.05;
    double upper_surface = .05*x*sqrt(4 - x*x);
    double camber_line   = .1*sin(M_PI*x/2.);
    double aux = SQ(upper_surface) - SQ(y-camber_line) + 0.002;
    if (inside_bounding_box && aux>0.0) 
      phi = pw_a*pow(aux,pw_n);
    return phi + potential_squarewell_yz(srch);
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
    else if (boost::iequals(efstr,"airfoil")) {
      potential = potential_airfoil;
      acceleration = acceleration_airfoil;
/*      for(double x= 0.0; x<2.0; x+=0.01) {
      for(double y=-0.5; y<0.5; y+=0.01) {

  {
    using namespace param;
    double phi = 0.0, ax = 0.0, ay = 0.0;
    assert (gdimension > 1);
    
    const double pw_n = extforce_sqwell_power;
    const double pw_a = extforce_sqwell_steepness;

    bool inside_bounding_box = std::abs(y)<0.25 && x>0 && x<2;
    double upper_surface = .05*x*sqrt(4 - x*x);
    double camber_line   = .1*sin(M_PI*x/2.);
    double aux = -SQ(y-camber_line) + SQ(upper_surface);
    if (inside_bounding_box && aux>0) 
      phi = pw_a*pow(aux,pw_n);
      ax = pw_n*pw_a*pow(phi,pw_n-1)
           * ( 2.*(y-camber_line)*(-.1*M_PI/2.*cos(M_PI/2.*x))
             - .0025*(4*x*(2 - x*x)));
      ay = pw_n*pw_a*pow(phi,pw_n-1) * 2.*(y-camber_line);
    std::cout << x << " " << y <<  " " << phi << " " << ax << " " << ay << " " << std::endl;
  }
      } std::cout << std::endl;
      }
  exit(0);
*/
    }
    else {
      clog_one(fatal) << "ERROR: bad external_force_type" << std::endl;
    }
  }

} // namespace external_force

#undef SQ
#endif // _eforce_h_
