/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
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
#include "density_profiles.h"
#define SQ(x) ((x)*(x))

namespace external_force {

  // acceleration and potential function types and pointers
  typedef double  (*potential_t)(const point_t&);
  typedef point_t (*acceleration_t)(const body&);
  static std::vector<potential_t> vec_potentials;
  static std::vector<acceleration_t> vec_accelerations;

  /**
   * @brief      1D walls: steep power-law-like potentials
   * @param      rp  Point coordinates
   */
  template <int I = 0>
  double potential_square_well(const point_t& rp) {
    using namespace param;
    const static double
       box[3] = {.5*box_length,.5*box_width,.5*box_height},
       pw_n = extforce_wall_powerindex,
       pw_a = extforce_wall_steepness;

    double phi = (((rp[I] < -box[I]) ? pow(-rp[I]- box[I], pw_n) : 0.0)
                 +((rp[I] >  box[I]) ? pow( rp[I]- box[I], pw_n) : 0.0))
                 *pw_a;
    return phi;
  }
  potential_t potential_walls_x = potential_square_well<0>;
  potential_t potential_walls_y = potential_square_well<1>;
  potential_t potential_walls_z = potential_square_well<2>;


  template <int I = 0>
  point_t acceleration_square_well(const body& particle) {
    using namespace param;
    point_t a = 0.0;
    point_t rp = particle.coordinates();
    const static double
       box[3] = {.5*box_length,.5*box_width,.5*box_height},
       pw_n = extforce_wall_powerindex,
       pw_a = extforce_wall_steepness;

    a[I]  = (((rp[I] <- box[I]) ? pow(-rp[I]- box[I], pw_n - 1) : 0.0)
            -((rp[I] >  box[I]) ? pow(rp[I] - box[I], pw_n - 1) : 0.0))
          *pw_n*pw_a;
    return a;
  }
  acceleration_t acceleration_walls_x = acceleration_square_well<0>;
  acceleration_t acceleration_walls_y = acceleration_square_well<1>;
  acceleration_t acceleration_walls_z = acceleration_square_well<2>;

  /**
   * @brief      Round or spherical boundary wall
   * @param      particle  The particle being accelerated
   */
  point_t acceleration_spherical_wall(const body& particle) {
    using namespace param;
    point_t a = 0.0;
    point_t rp = particle.coordinates();
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

  double potential_spherical_wall(const point_t& rp) {
    using namespace param;
    double phi = 0.0;
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
   * @brief      External force support for parabolic
   *             sphericall-symmetric density
   * @param      particle  The particle being accelerated
   */
  point_t acceleration_spherical_density_support (const body& particle) {
    using namespace param;
    point_t a = 0.0;
    static const double
        K0 = pressure_initial / pow(rho_initial, poly_gamma),
        rho0 = density_profiles::spherical_density_profile(0.);
    point_t rp = particle.coordinates();
    double r = rp[0]*rp[0];
    for (unsigned short i=1; i<gdimension; ++i)
      r += rp[i]*rp[i];
    r = sqrt(r);
    const double x = r / sphere_radius;
    if (x > 1e-12) {
      double rho = rho_initial / rho0
                 * density_profiles::spherical_density_profile(x);
      double drhodr = rho_initial / (rho0 * sphere_radius)
                    * density_profiles::spherical_drho_dr(x);
      double a_r = K0*poly_gamma*pow(rho,poly_gamma-2)*drhodr;
      for (short int i=0; i<gdimension; ++i)
        a[i] = a_r*rp[i] / r;
    }
    return a;
  }

  double potential_spherical_density_support(const point_t& rp) {
    using namespace param;
    static const double
        K0 = pressure_initial / pow(rho_initial, poly_gamma),
        rho0 = density_profiles::spherical_density_profile(0.);
    double r = rp[0]*rp[0];
    for (unsigned short i=1; i<gdimension; ++i)
      r += rp[i]*rp[i];
    r = sqrt(r);
    const double x = r / sphere_radius;
    double rho = rho_initial / rho0
               * density_profiles::spherical_density_profile(x);
    double phi = -K0*poly_gamma*pow(rho,poly_gamma-1.) / (poly_gamma-1.);
    return phi;
  }


  /**
   * @brief      Add uniform constant gravity acceleration
   * 	         in y-direction (or x-direction if number of
   * 	         dimensions == 1)
   * @param      particle  The particle being accelerated
   */
  point_t acceleration_gravity(const body& particle) {
    static const double grav = param::gravity_acceleration_constant;
    point_t acc = 0.0;
    if(gdimension > 1)
      acc[1] = -grav;  // negative y-direction
    else
      acc[0] = -grav;  // negative x-direction
    return acc;
  }

  double potential_gravity(const point_t& rp) {
    static const double grav = param::gravity_acceleration_constant;
    double height = rp[0];
    if(gdimension > 1)
      height = rp[1];
    return height*grav;
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
   * @param      particle  The particle being accelerated
   */
  point_t acceleration_airfoil (const body& particle) {
    using namespace param;
    point_t a = 0.0;
    assert (gdimension > 1);

    point_t rp =  particle.coordinates();
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
    return a;
  }

  double potential_airfoil (const point_t& rp) {
    using namespace param;
    double phi = 0.0;
    assert (gdimension > 1);

    static const double
                 alpha = airfoil_attack_angle*M_PI/180.0,
                 pw_n = extforce_wall_powerindex,
                 pw_a = extforce_wall_steepness;
    const double x1 = rp[0] - airfoil_anchor_x,
                 y1 = rp[1] - airfoil_anchor_y;
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
    return phi;
  }

  /**
   * @brief      Constant potential shift
   * @param      rp  Point coordinates
   */
  double potential_poison(const point_t& rp) {
    return param::zero_potential_poison_value;
  }

  /**
   * @brief      Total external force at a point 'srch'
   * @param      particle  Accelerated particle
   */
  point_t acceleration(const body& particle) {
    point_t a = 0.0;
    for (auto p : vec_accelerations)
      a += (*p)(particle);
    return a;
  }


  /**
   * @brief      Total external potential
   * @param      coords  Coordinates of where to compute the potential
   */
  double potential(const point_t& coords) {
    double phi = 0.0;
    for (auto p : vec_potentials)
      phi += (*p)(coords);
    return phi;
  }


  /**
   * @brief      External force selector
   * @param      efstr    ext. force string
   */
  void select(const std::string& efstr) {

    vec_potentials.clear();
    vec_accelerations.clear();

    if (boost::iequals(efstr,"zero") or boost::iequals(efstr,"none"))
      return; // trivial case

    using namespace std;

    // parse efstr: external force specification string is a comma-separated
    // list of potentials / accelerations which need to be added up: e.g.
    // "spherical wall,walls:xyz,gravity"
    vector<string> split_efstr;
    boost::split(split_efstr, efstr, boost::is_any_of(","));
    for (auto it = split_efstr.begin(); it!= split_efstr.end(); ++it) {
      if (boost::iequals(*it,"spherical wall")) {
        vec_potentials.push_back(potential_spherical_wall);
        vec_accelerations.push_back(acceleration_spherical_wall);
      }
      else if (boost::iequals(*it,"airfoil")) {
        vec_potentials.push_back(potential_airfoil);
        vec_accelerations.push_back(acceleration_airfoil);
      }
      else if (boost::iequals(*it,"spherical density support")) {
        density_profiles::select();
        vec_potentials.push_back(potential_spherical_density_support);
        vec_accelerations.push_back(acceleration_spherical_density_support);
      }
      else if (boost::iequals(*it,"gravity")) {
        vec_potentials.push_back(potential_gravity);
        vec_accelerations.push_back(acceleration_gravity);
      }
      else if (boost::iequals(it->substr(0,6),"walls:")) {
        // parse in which directions to place the walls
        // this can be e.g. "walls:xyz" or "walls:y" etc.
        const char *cxyz = it->substr(6).c_str();
        char   imx = min(3,(int)it->substr(6).length());
        for (int i=0; i<imx; ++i) {
          switch (cxyz[i]) {
          case 'x':
          case 'X':
            vec_potentials.push_back(potential_walls_x);
            vec_accelerations.push_back(acceleration_walls_x);
            break;
          case 'y':
          case 'Y':
            vec_potentials.push_back(potential_walls_y);
            vec_accelerations.push_back(acceleration_walls_y);
            break;
          case 'z':
          case 'Z':
            vec_potentials.push_back(potential_walls_z);
            vec_accelerations.push_back(acceleration_walls_z);
            break;
          default:
            clog_fatal("ERROR: bad external_force_type" << std::endl);
            assert(false);
          }
        }
      }
      else if (boost::iequals(*it,"poison")) {
        // zero potential shift
        vec_potentials.push_back(potential_poison);
      }
      else {
        clog_fatal("ERROR: bad external_force_type" << std::endl);
      }
    } // for it in split_efstr

  } // select()


  /**
   * @brief      Artificial drag force - used for
   *             particle relaxation
   * @param      vel   Velocity against the drag
   */
  point_t acceleration_drag(const point_t& vel) {
    using namespace param;
    point_t acc = 0.0;
    double v2 = vel[0]*vel[0];
    for (short int i=1; i<gdimension; ++i)
      v2 += vel[i]*vel[i];

    acc -= (relaxation_beta + relaxation_gamma*v2) * vel;
    return acc;
  }


} // namespace external_force

#undef SQ
#endif // _eforce_h_
