/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file wvt.h
 * @author I. Sagert
 * @date August 2019
 * @brief WVT relaxation 
 */

#ifndef _wvt_h_
#define _wvt_h_

#include <vector>

#include "params.h"
#include "utils.h"
#include "user.h"
#include "kernels.h"
#include "tree.h"

#include "eos.h"
#include "integration.h"


namespace wvt{
  using namespace param;

  /** 
   * @brief:    Converts cartesian particle coordinates into 
   *            spherical 
   *
   * @param     pos_c   Particle position in cartesian coord.
   */
  point_t
  cartesian_to_spherical (const point_t& pos_c) {
    double r = norm2(pos_c);
    point_t pos_s = 0.0;

    if (gdimension == 1) {
      pos_s[0] = r;
    }
    else if (gdimension == 2) {
      pos_s[0] = r;
      pos_s[1] = atan2(pos_c[1],pos_c[0]);
    }
    else {
      pos_s[0] = r;
      pos_s[1] = atan2(pos_c[1],pos_c[0]);
      pos_s[2] = acos(pos_c[2]/r);
    }
    return pos_s;
  }


  /** 
   * @brief:    Converts spherical particle coordinates into 
   *            cartesian 
   *
   * @param     pos_s   Particle position in spherical coord.
   */
  point_t
  spherical_to_cartesian (const point_t& pos_s) {
    double r = pos_s[0];
    point_t pos_c = 0.0;

    if (gdimension == 1) {
      pos_c[0] = r;
    }
    else if (gdimension == 2) {
      pos_c[0] = r*cos(pos_s[1]);
      pos_c[1] = r*sin(pos_s[1]);
    }
    else {
      pos_c[0] = r*sin(pos_s[2])*cos(pos_s[1]);
      pos_c[1] = r*sin(pos_s[2])*sin(pos_s[1]);
      pos_c[2] = r*cos(pos_s[2]);
    }
    return pos_c;
  }


  /**
   * @brief      Pseudo-acceleration for WVT equilibration 
   *             [Arth et al. 2019]
   *
   * @param      particle  The particle body
   * @param      nbs       Vector of neighbor particles
   */
  void
  wvt_acceleration_arth(
      body& particle,
      std::vector<body*>& nbs)
  {
    using namespace param;
    using namespace kernels;

    particle.setNeighbors(nbs.size());

    // this particle (index 'a')
    const point_t pos_a = particle.coordinates();
    const double r_a    = norm2(pos_a);
    const double h_a    = particle.radius();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();

    point_t acc_a = 0.0;

    // Set particle acceleration only inside the sphere
    if (r_a <= sphere_radius) {
      
      // Loop over all neighbors and calculate repulsive forces
      for(int b = 0; b < n_nb; ++b) {
        const body * const nb = nbs[b];
        double h_b  = nb->radius();
        double h_ab = 0.5*(h_a + h_b);

        point_t pos_b = nb->coordinates();
        double r_ab = flecsi::distance(pos_a,pos_b);
        double W_ab = sph_kernel_function(r_ab,h_ab)*pow(h_ab,gdimension);
        if (r_ab > 0.0 && r_ab <= h_ab) {
          acc_a += (h_ab*W_ab/r_ab) * (pos_a-pos_b);
        }

        // If particle sees the edge, interact with mirror particles 
        if ((r_a + h_a) > sphere_radius) {
          point_t pos_bs = cartesian_to_spherical(pos_b);
          if (pos_bs[0] < sphere_radius) {
            pos_bs[0] = sphere_radius + (sphere_radius - pos_bs[0]);
            pos_b = spherical_to_cartesian(pos_bs);
            r_ab  = flecsi::distance(pos_a,pos_b);
            W_ab  = sph_kernel_function(r_ab,h_ab)*pow(h_ab,gdimension);
            if (r_ab > 0.0 && r_ab <= h_ab) {
              acc_a += (h_ab*W_ab/r_ab) * (pos_a-pos_b);
            }
          }
        }
      }

      // If particle sees the edge, interact with own mirror image 
      if ((r_a + h_a) > sphere_radius) {
        point_t pos_ms = cartesian_to_spherical(pos_a);
        pos_ms[0] = sphere_radius + (sphere_radius - pos_ms[0]);

        point_t pos_m = spherical_to_cartesian(pos_ms);
        double r_am = flecsi::distance(pos_a,pos_m);
        double W_am = sph_kernel_function(r_am,h_a)*pow(h_a,gdimension);
        if (r_am > 0.0 && r_am <= h_a) {
          acc_a += (h_a*W_am/r_am) * (pos_a-pos_m);
        }
      }
    }
    particle.setAcceleration(acc_a);
  }//wvt_acceleration_arth



  /**
   * @brief      Pseudo-acceleration for WVT equilibration 
   *             [Diehl, PASA 2015]
   *
   * @param      particle  The particle body
   * @param      nbs       Vector of neighbor particles
   */
  void
  wvt_acceleration_diehl(
      body& particle,
      std::vector<body*>& nbs)
  {
    using namespace param;
    using namespace kernels;

    particle.setNeighbors(nbs.size());

    // this particle (index 'a')
    const point_t pos_a = particle.coordinates();
    const double r_a    = norm2(pos_a);
    const double h_a    = particle.radius();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();

    point_t acc_a = 0.0;

    // Set particle acceleration only inside the sphere
    if (r_a <= sphere_radius) {

      // Loop over all neighbors and calculate repulsive forces
      for(int b = 0; b < n_nb; ++b) {
        const body * const nb = nbs[b];
        double h_b  = nb->radius();
        double h_ab = 0.5*(h_a + h_b);

        point_t pos_b = nb->coordinates();
        double r_ab = flecsi::distance(pos_a,pos_b);
        double W_ab = h_ab*h_ab/((r_ab+0.3*h_ab)*(r_ab+0.3*h_ab)) 
                    - h_ab*h_ab/((h_ab+0.3*h_ab)*(h_ab+0.3*h_ab));
        if (r_ab > 0.0 && r_ab <= h_ab) {
           acc_a += (h_ab*W_ab) * (pos_a-pos_b)/r_ab;
        }
        
        // If particle sees the edge, interact with mirror particles 
        if ((r_a + h_a) > sphere_radius) {
          point_t pos_bs = cartesian_to_spherical(pos_b);
          if (pos_bs[0] < sphere_radius) {
            pos_bs[0] = sphere_radius + (sphere_radius - pos_bs[0]);
            pos_b = spherical_to_cartesian(pos_bs);
            r_ab  = flecsi::distance(pos_a,pos_b);
            W_ab  = h_ab*h_ab/((r_ab+0.3*h_ab)*(r_ab+0.3*h_ab)) 
                  - h_ab*h_ab/((h_ab+0.3*h_ab)*(h_ab+0.3*h_ab));
            if (r_ab > 0.0 && r_ab <= h_ab) {
              acc_a += (h_ab*W_ab) * (pos_a-pos_b)/r_ab;
            }
          }
        }
      }

      // If particle sees the edge, interact with own mirror image 
      if ((r_a + h_a) > sphere_radius) {
        point_t pos_ms = cartesian_to_spherical(pos_a);
        pos_ms[0] = sphere_radius + (sphere_radius - pos_ms[0]);

        point_t pos_m = spherical_to_cartesian(pos_ms);
        double r_am = flecsi::distance(pos_a,pos_m);
        double W_am = h_a*h_a/((r_am+0.3*h_a)*(r_am+0.3*h_a)) 
                    - h_a*h_a/((h_a+0.3*h_a)*(h_a+0.3*h_a));
        if (r_am > 0.0 && r_am <= h_a) {
          acc_a += (h_a*W_am) * (pos_a-pos_m)/r_am;
        }
      }
    }
    particle.setAcceleration(acc_a);
  }//wvt_acceleration_diehl



  /**
   * @brief: Rescales smoothing length following prescription of 
   *         Diehl et al., PASA 2015. Does not seem to work yet. 
   */
  void
  compute_smoothinglength_wvt(
      std::vector<body>& bodies)
  {
    double Vsph = 0.0;
    if (gdimension == 1) {
      for(size_t i = 0; i < bodies.size(); ++i){
        double h_a = bodies[i].radius();
        Vsph += 2.0*h_a;
      }
      double Vtotal  = 2.0*sphere_radius;
      double scaling = Vtotal*wvt_ngb/Vsph;
      #pragma omp parallel for
      for(size_t i = 0; i < bodies.size(); ++i){
        double h_a = bodies[i].radius();
        bodies[i].set_radius(scaling * h_a);
      }
    }
    else if (gdimension == 2) {
      for(size_t i = 0 ; i < bodies.size(); ++i){
        double h_a = bodies[i].radius();
        Vsph += M_PI*h_a*h_a;
      }
      double Vtotal  = M_PI*sphere_radius*sphere_radius;
      double scaling = sqrt(Vtotal*wvt_ngb/Vsph);
      #pragma omp parallel for
      for(size_t i = 0 ; i < bodies.size(); ++i){
        double scaling = sqrt(Vtotal*wvt_ngb/Vsph);
        double h_a = bodies[i].radius();
        bodies[i].set_radius(scaling * h_a);
      }
    }
    else {
      for(size_t i = 0 ; i < bodies.size(); ++i){
        double h_a = bodies[i].radius();
        Vsph += M_PI*4.0*h_a*h_a*h_a/3.0;
      }
      double Vtotal = 4.0*M_PI*sphere_radius*sphere_radius
                    * sphere_radius/3.0;
      double scaling = pow(Vtotal*wvt_ngb/Vsph,1.0/3.0);
      #pragma omp parallel for
      for(size_t i = 0 ; i < bodies.size(); ++i){
        double h_a = bodies[i].radius()*scaling;
        bodies[i].set_radius(h_a);
      }
    }
  } //compute_smoothinglength_wvt


  /**
   * @brief      WVT displacement  dx
   *             [Diehl et al., PASA 2015]
   *
   *             dx = mu*h_a*sum_b[f(h_ab,r_ab)*r_ab]
   *
   * @param      srch  The source's body holder
   */
  void
  wvt_displacement (body& source) {
    using namespace param;    
    double wvt_mu_it = wvt_mu;

    // Check where particle will end up in next iteration
    point_t rp = source.coordinates() + wvt_mu_it
               * source.getAcceleration();

    double mass = source.mass();
    double r = norm2(rp);

    // Always decrease "timestep" when particles move too far out.
    // This should be global. Surprisingly it also works when the 
    // stepsize is updated only locally. 
    while (r/sphere_radius > 1.2) {
        wvt_mu_it *= 0.5;
        rp = source.coordinates() + wvt_mu_it
           * source.getAcceleration();
        r = norm2(rp);
    }

    // Freeze particels in the outer edge 
    if (boost::iequals(wvt_boundary, "frozen")) {
      point_t pos = source.coordinates();
      r = norm2(pos);
      if (r/sphere_radius >= 0.9) {
        source.setAcceleration(0.0);
        rp = source.coordinates();
      }
    }
    // Reflect particles on the sphere edge 
    else if (boost::iequals(wvt_boundary, "reflective")) {
      if (r/sphere_radius >= 1.0) {
        point_t rs = cartesian_to_spherical(rp);
        rs[0] = rs[0] - 2.0*(rs[0] - sphere_radius);
        rp = spherical_to_cartesian(rs);
      }
    }
    else {
      std::cout << "Error: wvt boundary undefined!" << std::endl;
      exit(0);
    }

    r = norm2(rp);
    double rho0 = density_profiles::spherical_density_profile(0);
    double rho  = ((param::rho_initial)/rho0) 
            * density_profiles::spherical_density_profile(r/sphere_radius);
    double h_a = kernels::kernel_width*pow(mass/rho,1./gdimension);
    source.set_radius(h_a);
    source.set_coordinates(rp);
    source.setDensity(rho);
  } // wvt_displacement


  // wvt types and pointers
  typedef void (*compute_quantity_t)(body&, std::vector<body*>&);
  compute_quantity_t wvt_acceleration = wvt_acceleration_diehl;

  /**
   * @brief      WVT method selector
   */
  void select() {
    using namespace param;
    if (boost::iequals(wvt_method,"diehl")) {
      wvt_acceleration = wvt_acceleration_diehl;
    }
    else if (boost::iequals(wvt_method,"arth")) {
      wvt_acceleration = wvt_acceleration_arth;
    }
    else {
      clog(error) << "ERROR: wrong parameter in wvt";
      exit(2);
    }
  } // select()


}; //wvt


 #endif // _wvt_h_
