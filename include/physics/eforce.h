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
    const double pw_n = extforce_wall_powerindex;
    const double pw_a = extforce_wall_steepness;
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
    const double pw_n = extforce_wall_powerindex;
    const double pw_a = extforce_wall_steepness;
    for (unsigned short i=1; i<gdimension; ++i) {
      phi += (((rp[i] <- box[i]) ? pow(-rp[i]- box[i], pw_n) : 0.0)
             +((rp[i] >  box[i]) ? pow(rp[i] - box[i], pw_n) : 0.0))
            *pw_a;
    }
    return phi;
  }


  /**
   * @brief      Round or spherical boundary wall
   * @param      srch  The source's body holder
   */
  point_t acceleration_spherical_wall(body_holder* srch) {
    using namespace param;
    point_t a = 0.0;
    point_t rp =  srch->getBody()->getPosition();
    const double pw_n = extforce_wall_powerindex;
    const double pw_a = extforce_wall_steepness;
    double r = rp[0]*rp[0];
    for (unsigned short i=1; i<gdimension; ++i) 
      r += rp[i]*rp[i];
    r = sqrt(r);
    if (r > sphere_radius) {
      const double ar = pw_n*pw_a*pow(r - sphere_radius, pw_n - 1);
      for (unsigned short i=0; i<gdimension; ++i) 
        a[i] = -rp[i]/r * ar;
    }
     
    return a;
  }

  double potential_spherical_wall(body_holder* srch) {
    using namespace param;
    double phi = 0.0;
    point_t rp =  srch->getBody()->getPosition();
    const double pw_n = extforce_wall_powerindex;
    const double pw_a = extforce_wall_steepness;
    double r = rp[0]*rp[0];
    for (unsigned short i=1; i<gdimension; ++i) 
      r += rp[i]*rp[i];
    r = sqrt(r);
    if (r > sphere_radius) 
      phi = pw_a*pow(r - sphere_radius, pw_n);
    return phi;
  }


  /**
   * @brief      2D airfoil in a wind tunnel
   *
   * The airfoil profile is centered at the anchor, tilted
   * at an angle to the flow. The shape of the airfoil can be described by the
   * following three parameters:
   *  - airfoil_size:           airfoil horizontal extent;
   *  - airfoil_thickness:      how thick is it;
   *  - airfoil_camber:         maximum deviation of camber line from the chord.
   * 
   * Airfoil is positioned and rotated relative to its rear tip:
   *  - airfoil_anchor_x:       the x-coordinate of the anchor;
   *  - airfoil_anchor_y:       the y-coordinate of the anchor;
   *  - airfoil_attack_angle:   angle of attack - rotation from initial position
   *                            which is parallel to the x-axis.
   *
   * @param      srch  The source's body holder
   */
  point_t acceleration_airfoil (body_holder* srch) {
    using namespace param;
    point_t a = 0.0;
    assert (gdimension > 1);
    
    point_t rp =  srch->getBody()->getPosition();
    const double x1 = rp[0] - airfoil_anchor_x,
                 y1 = rp[1] - airfoil_anchor_y,
                 alpha = airfoil_attack_angle*M_PI/180.0,
                 pw_n = extforce_wall_powerindex,
                 pw_a = extforce_wall_steepness;
    const double x = x1*cos(alpha) + y1*sin(alpha),
                 y =-x1*sin(alpha) + y1*cos(alpha);

    bool inside_bounding_box = std::abs(y)<5.0*airfoil_thickness
        && x>-airfoil_size*0.02 && x< airfoil_size*1.02;
    double upper_surface = airfoil_thickness*x
                         * sqrt(airfoil_size*airfoil_size - x*x);
    double camber_line   = airfoil_camber*sin(M_PI*x/2.);
    double phi = SQ(upper_surface) - SQ(y-camber_line) + 0.002;
    if (inside_bounding_box && phi>0.0) {
      double a0, a1;
      a0 = pw_n*pw_a*pow(phi,pw_n-1)
           * ( 2.*(y-camber_line)*(-airfoil_camber*M_PI/2.*cos(M_PI/2.*x))
             - airfoil_thickness*airfoil_thickness*2*x
                  *(airfoil_size*airfoil_size  -   2*x*x));
      a1 = pw_n*pw_a*pow(phi,pw_n-1) * 2.*(y-camber_line);
      a[0] = a0*cos(alpha) - a1*sin(alpha);
      a[1] = a0*sin(alpha) + a1*cos(alpha);
    }
    return a + acceleration_squarewell_yz(srch);
  }

  double potential_airfoil (body_holder* srch)
  {
    using namespace param;
    double phi = 0.0;
    assert (gdimension > 1);
    
    point_t rp =  srch->getBody()->getPosition();
    const double x1 = rp[0] - airfoil_anchor_x,
                 y1 = rp[1] - airfoil_anchor_y,
                 alpha = airfoil_attack_angle*M_PI/180.0,
                 pw_n = extforce_wall_powerindex,
                 pw_a = extforce_wall_steepness;
    const double x = x1*cos(alpha) + y1*sin(alpha),
                 y =-x1*sin(alpha) + y1*cos(alpha);

    bool inside_bounding_box = std::abs(y)<5.0*airfoil_thickness
        && x>-airfoil_size*0.02 && x< airfoil_size*1.02;
    double upper_surface = airfoil_thickness*x
                         * sqrt(airfoil_size*airfoil_size - x*x);
    double camber_line   = airfoil_camber*sin(M_PI*x/2.);
    double aux = SQ(upper_surface) - SQ(y-camber_line) + 0.002;
    if (inside_bounding_box && aux>0.0) 
      phi = pw_a*pow(aux,pw_n);
    return phi + potential_squarewell_yz(srch);
  }

  /**
   * @brief      Drag case : apply drag force to 
   * 	         relax the initial star during initial
   * 	         few steps
   * @param      srch  The source's body holder
   */
  point_t acceleration_do_drag(body_holder* srch) {
    using namespace param;
    body* source = srch->getBody();
    int64_t iteration = 0;
    point_t a = 0.0;
    if(do_drag && iteration <= relax_steps){
      //Redefine drag coefficient with dt
      double drag_coeff_dt = drag_coeff/initial_dt;
      a -=drag_coeff_dt*source->getVelocity();
    }
    return a;
  }

  double potential_do_drag(body_holder* srch) {
    // No potential for dragging. Only du/dt = 0 during relaxation
    return 0.0; // POISON IT
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
    else if (boost::iequals(efstr,"spherical wall")) {
      potential = potential_spherical_wall;
      acceleration = acceleration_spherical_wall;
    }
    else if (boost::iequals(efstr,"airfoil")) {
      potential = potential_airfoil;
      acceleration = acceleration_airfoil;
    }
    else if (boost::iequals(efstr,"drag")) {
      potential = potential_do_drag;
      acceleration = acceleration_do_drag;
    }
    else {
      clog_one(fatal) << "ERROR: bad external_force_type" << std::endl;
    }
  }

} // namespace external_force

#undef SQ
#endif // _eforce_h_
