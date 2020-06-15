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

#include "params.h"
#include <vector>

// Basic elements for the wvt
namespace wvt_basic {
using namespace param;
bool wvt_converged = false;
} // namespace wvt_basic

#include "kernels.h"
#include "tree.h"
#include "user.h"
#include "utils.h"

#include "eos.h"
#include "integration.h"

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))

namespace wvt {
using namespace param;

//  /**
//   * @brief: Resets WVT sphere boundary for
//   *         based on h, for better stability
//   */
//  void
//  wvt_set_radius(
//      std::vector<body>& bodies)
//  {
//    int n_r = 1000;
//    double dr = sphere_radius/(n_r*1.0);
//    double mass = 0.0;
//
//    for(auto& b: bodies) {
//      mass = b.mass();
//    }
//
//    double h = 0.0;
//    double r = sphere_radius;
//    int i = 0;
//
//    while(i<=n_r && h < wvt_boundary_scale*sphere_radius) {
//      i++;
//      r = dr*i;
//      double rho0 = density_profiles::spherical_density_profile(0);
//      double rho  = (param::rho_initial)/rho0
//                  * density_profiles::spherical_density_profile(r/sphere_radius);
//      h = sph_eta*kernels::kernel_width*pow(mass/rho,1./gdimension);
//    }
//    wvt_basic::wvt_radius = r;
//  } // reset_wvt_boundary

void
compute_density(body & particle, std::vector<body *> & nbs) {
  using namespace kernels;
  const double h_a = particle.radius();
  const point_t pos_a = particle.coordinates();
  const int n_nb = nbs.size();
  mpi_assert(n_nb > 0);

  double rho_a = 0.0;
  for(int b = 0; b < n_nb; ++b) { // Vectorized
    const body * const nb = nbs[b];
    double m_b = nb->mass();
    double h_b = nb->radius();
    point_t pos_b = nb->coordinates();
    double r_ab = flecsi::magnitude(pos_a - pos_b);
    double Wab = sph_kernel_function(r_ab, .5 * (h_a + h_b));
    rho_a += m_b * Wab;
  } // for

  if(not(rho_a > 0)) {
    std::cout << "Density of a particle is not a positive number: "
              << "rho = " << rho_a << std::endl;
    std::cout << "Failed particle id: " << particle.id() << std::endl;
    std::cerr << "particle position: " << particle.coordinates() << std::endl;
    std::cerr << "particle velocity: " << particle.getVelocity() << std::endl;
    std::cerr << "particle acceleration: " << particle.getAcceleration()
              << std::endl;
    std::cerr << "smoothing length:  " << particle.radius() << std::endl;
    assert(false);
  }
  particle.setDensity(rho_a);
}

/**
 * @brief:    Converts cartesian particle coordinates into
 *            spherical
 *
 * @param     pos_c   Particle position in cartesian coord.
 */
point_t
cartesian_to_spherical(const point_t & pos_c) {
  double r = magnitude(pos_c);
  point_t pos_s = 0.0;

  if(gdimension == 1) {
    pos_s[0] = r;
  }
  else if(gdimension == 2) {
    pos_s[0] = r;
    pos_s[1] = atan2(pos_c[1], pos_c[0]);
  }
  else {
    pos_s[0] = r;
    pos_s[1] = atan2(pos_c[1], pos_c[0]);
    pos_s[2] = acos(pos_c[2] / r);
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
spherical_to_cartesian(const point_t & pos_s) {
  double r = pos_s[0];
  point_t pos_c = 0.0;

  if(gdimension == 1) {
    pos_c[0] = r;
  }
  else if(gdimension == 2) {
    pos_c[0] = r * cos(pos_s[1]);
    pos_c[1] = r * sin(pos_s[1]);
  }
  else {
    pos_c[0] = r * sin(pos_s[2]) * cos(pos_s[1]);
    pos_c[1] = r * sin(pos_s[2]) * sin(pos_s[1]);
    pos_c[2] = r * cos(pos_s[2]);
  }
  return pos_c;
}

/**
 * @brief:    Removes particle acceleration component that
 *            is parallel to the direction of motion as a
 *            way to implement frozen boundary conditions.
 *
 * @param     pos_s   Particle position in cartesian coord.
 *            acc_c   Particle acceleration
 */
point_t
remove_radial_acc(const point_t & pos_c, const point_t & acc_c) {
  double r = magnitude(pos_c);
  point_t pos_n = (1.0 / r) * pos_c;
  point_t acc_n = acc_c - dot(pos_n, acc_c) * pos_n;
  return 2.0 * acc_n;
}

/**
 * @brief      Pseudo-acceleration for WVT equilibration
 *             [Arth et al. 2019]
 *
 * @param      particle  The particle body
 * @param      nbs       Vector of neighbor particles
 */
void
wvt_acceleration_arth(body & particle, std::vector<body *> & nbs) {
  using namespace param;
  using namespace kernels;

  particle.setNeighbors(nbs.size());

  // this particle (index 'a')
  const point_t pos_a = particle.coordinates();
  double r_a = magnitude(pos_a);
  double h_a = particle.radius();

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  point_t acc_a = 0.0;
  double boundary = wvt_radius;

  // Set particle acceleration only inside the sphere
  if(r_a <= boundary) {
    // Loop over all neighbors and calculate repulsive forces
    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      double h_b = nb->radius();
      point_t pos_b = nb->coordinates();
      double r_ab = flecsi::distance(pos_a, pos_b);
      double h_ab = 0.5 * (h_a + h_b);
      double W_ab = sph_kernel_function(r_ab, h_ab) * pow(h_ab, gdimension);
      // if (r_ab > 0.0 && r_ab <= h_ab) {
      if(r_ab > 0.0) {
        acc_a += (h_ab * W_ab) * (pos_a - pos_b) / r_ab;
      }
      if(boost::iequals(wvt_boundary, "reflective")) {
        // If particle sees the edge, interact with mirror particles
        if((r_a + h_a) > boundary) {
          point_t pos_bs = cartesian_to_spherical(pos_b);
          if(pos_bs[0] < boundary) {
            pos_bs[0] = boundary + (boundary - pos_bs[0]);
            pos_b = spherical_to_cartesian(pos_bs);
            r_ab = flecsi::distance(pos_a, pos_b);
            // W_ab  = sph_kernel_function(r_ab,h_ab)*pow(h_ab,gdimension);
            // if (r_ab > 0.0 && r_ab <= h_ab) {
            //  acc_a += (h_ab*W_ab)*(pos_a-pos_b)/r_ab;
            W_ab = sph_kernel_function(r_ab, h_ab) * pow(h_a, gdimension);
            // if (r_ab > 0.0 && r_ab <= h_a) {
            if(r_ab > 0.0) {
              acc_a += (h_a * W_ab) * (pos_a - pos_b) / r_ab;
            }
          }
        }
      }
    } // Loop over neighbors

    // If particle sees the edge, interact with own mirror image
    if(boost::iequals(wvt_boundary, "reflective")) {
      if((r_a + h_a) > boundary) {
        point_t pos_ms = cartesian_to_spherical(pos_a);
        pos_ms[0] = boundary + (boundary - pos_ms[0]);
        point_t pos_m = spherical_to_cartesian(pos_ms);
        double r_am = flecsi::distance(pos_a, pos_m);
        double W_am = sph_kernel_function(r_am, h_a) * pow(h_a, gdimension);
        // if (r_am > 0.0 && r_am <= h_a) {
        if(r_am > 0.0) {
          acc_a += (h_a * W_am) * (pos_a - pos_m) / r_am;
        }
      }
    } // reflective bc check
    // Remove acceleration in radial direction if particle sees edge
    if(boost::iequals(wvt_boundary, "frozen")) {
      if((r_a + h_a) > boundary) {
        acc_a = remove_radial_acc(pos_a, acc_a);
      }
    } // frozen bc check
  }
  particle.setAcceleration(acc_a);
} // wvt_acceleration_arth

/**
 * @brief      Pseudo-acceleration for WVT equilibration
 *             following Diehl et al., PASA 2015
 *
 * @param      particle  The particle body
 * @param      nbs       Vector of neighbor particles
 */
void
wvt_acceleration_diehl(body & particle, std::vector<body *> & nbs) {
  using namespace param;
  using namespace kernels;

  particle.setNeighbors(nbs.size());

  // this particle (index 'a')
  const point_t pos_a = particle.coordinates();
  const double r_a = magnitude(pos_a);
  const double h_a = particle.radius();
  const double eps = 0.3;

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  point_t acc_a = 0.0;
  double boundary = wvt_radius;

  // Set particle acceleration only inside the sphere
  if(r_a <= boundary) {
    // Loop over all neighbors and calculate repulsive forces
    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      double h_b = nb->radius();
      double h_ab = 0.5 * (h_a + h_b);

      point_t pos_b = nb->coordinates();
      double r_ab = flecsi::distance(pos_a, pos_b);
      double W_ab =
        SQ(h_ab / (r_ab + eps * h_ab)) - SQ(h_ab / (h_ab + eps * h_ab));
      if(r_ab > 0.0) {
        double W_abt = std::max(W_ab, 0.0);
        acc_a += (h_a * W_abt) * (pos_a - pos_b) / r_ab;
      }
      if(boost::iequals(wvt_boundary, "reflective")) {
        // If particle sees the edge, interact with mirror particles
        if((r_a + h_a) > boundary) {
          point_t pos_bs = cartesian_to_spherical(pos_b);
          if(pos_bs[0] < boundary) {
            pos_bs[0] = boundary + (boundary - pos_bs[0]);
            pos_b = spherical_to_cartesian(pos_bs);
            r_ab = flecsi::distance(pos_a, pos_b);
            // W_ab  = SQ(h_ab/(r_ab+eps*h_ab))-SQ(h_ab/(h_ab+eps*h_ab));
            W_ab = SQ(h_a / (r_ab + eps * h_a)) - SQ(h_a / (h_a + eps * h_a));
            if(r_ab > 0.0) {
              double W_abt = std::max(W_ab, 0.0);
              // acc_a += (h_ab*W_abt)*(pos_a-pos_b)/r_ab;
              acc_a += (h_a * W_abt) * (pos_a - pos_b) / r_ab;
            }
          }
        }
      } // reflective bc check
    } // Loop over neighbors
    if(boost::iequals(wvt_boundary, "reflective")) {
      // If particle sees the edge, interact with own mirror image
      if((r_a + h_a) > boundary) {
        point_t pos_ms = cartesian_to_spherical(pos_a);
        pos_ms[0] = boundary + (boundary - pos_ms[0]);

        point_t pos_m = spherical_to_cartesian(pos_ms);
        double r_am = flecsi::distance(pos_a, pos_m);
        double W_am =
          SQ(h_a / (r_am + eps * h_a)) - SQ(h_a / (h_a + eps * h_a));
        if(r_am > 0.0) {
          double W_amt = std::max(W_am, 0.0);
          acc_a += (h_a * W_amt) * (pos_a - pos_m) / r_am;
        }
      }
    } // reflective bc check
    if(boost::iequals(wvt_boundary, "frozen")) {
      if((r_a + h_a) > boundary) {
        acc_a = remove_radial_acc(pos_a, acc_a);
      }
    } // frozen bc check
  }
  particle.setAcceleration(acc_a);
} // wvt_acceleration_diehl

/**
 * @brief: Rescales smoothing length according to
 *         total SPH particle count and desirted
 *         number of neighbors
 *         See e.g. Diehl et al., PASA 2015.
 */
void
compute_smoothinglength_wvt(std::vector<body> & bodies) {
  double Vsph = 0.0;
  double scaling = 1.0;
  double boundary = wvt_radius;
  if(gdimension == 1) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      double h = bodies[i].radius();
      Vsph += 2.0 * h;
    }
    mpi_utils::reduce_sum(Vsph);
    double Vtotal = 2.0 * boundary;
    scaling = Vtotal * wvt_ngb / Vsph;
  }
  else if(gdimension == 2) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      double h = bodies[i].radius();
      Vsph += M_PI * SQ(h);
    }
    mpi_utils::reduce_sum(Vsph);
    double Vtotal = M_PI * SQ(boundary);
    scaling = sqrt(Vtotal * wvt_ngb / Vsph);
  }
  else {
    for(size_t i = 0; i < bodies.size(); ++i) {
      double h = bodies[i].radius();
      Vsph += M_PI * 4.0 * CU(h) / 3.0;
    }
    mpi_utils::reduce_sum(Vsph);
    double Vtotal = 4.0 * M_PI * CU(boundary) / 3.0;
    scaling = pow(Vtotal * wvt_ngb / Vsph, 1.0 / 3.0);
  }
  for(size_t i = 0; i < bodies.size(); ++i) {
    point_t pos = bodies[i].coordinates();
    double r_a = magnitude(pos);
    double h = bodies[i].radius() * scaling;
    bodies[i].set_radius(h);
  }
} // compute_smoothinglength_wvt

/**
 * @brief: Checks convergence of WVT
 *         Stops the relaxation as soon as the most
 *         particles are moved less than a small fraction
 *         of inter-particle distance.
 *         See Arth et al. 2019 for details
 */
void
check_convergence_wvt(std::vector<body> & bodies) {
  int cnt = 0, cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;
  int nbs = bodies.size();

  for(auto & b : bodies) {
    // distance that particles moved
    point_t delta = wvt_mu * b.getAcceleration();
    const double d = magnitude(delta);
    const double h = b.radius();

    // inter-particle distance, depending on dimension
    double d_mps = pow(1.0 / wvt_ngb, 1.0 / gdimension) * h;

    //      if (boost::iequals(wvt_method,"arth")) {
    //        if (gdimension == 1) {
    //          d_mps = h/wvt_ngb;
    //        }
    //        else if (gdimension == 2) {
    //          d_mps = sqrt(M_PI/wvt_ngb)*h;
    //        }
    //        else {
    //          d_mps = pow((4.0*M_PI/3.0)/wvt_ngb,1.0/3.0)*h;
    //        }
    //      }

    if(d > 1.0 * d_mps) {
      ++cnt;
    }
    if(d > 0.1 * d_mps) {
      ++cnt1;
    }
    if(d > 0.01 * d_mps) {
      ++cnt2;
    }
    if(d > 0.001 * d_mps) {
      ++cnt3;
    }
    if(d > 0.0001 * d_mps) {
      ++cnt4;
    }
  } // Loop over bodies

  mpi_utils::reduce_sum(nbs);
  mpi_utils::reduce_sum(cnt);
  mpi_utils::reduce_sum(cnt1);
  mpi_utils::reduce_sum(cnt2);
  mpi_utils::reduce_sum(cnt3);
  mpi_utils::reduce_sum(cnt4);

  double moveMPS[5];
  moveMPS[0] = cnt * 100.0 / (nbs * 1.0);
  moveMPS[1] = cnt1 * 100.0 / (nbs * 1.0);
  moveMPS[2] = cnt2 * 100.0 / (nbs * 1.0);
  moveMPS[3] = cnt3 * 100.0 / (nbs * 1.0);
  moveMPS[4] = cnt4 * 100.0 / (nbs * 1.0);

  log_one(trace) << "consider converged when " << moveMPS[3] << " < "
                 << wvt_convergence_point << std::endl;
  if(moveMPS[3] < wvt_convergence_point) {
    wvt_basic::wvt_converged = true;
  }
} // check convergence

/**
 * @brief: Calculates standard deviation of particle
 *         densities from target density
 */
void
calculate_standard_deviation(std::vector<body> & bodies) {
  double nominator = 0.0;
  int nbs = bodies.size();

  for(auto & b : bodies) {
    point_t pos = b.coordinates();
    double r = magnitude(pos);
    double rho0 = density_profiles::spherical_density_profile(0);
    double rho_wd =
      density_profiles::spherical_density_profile(r / sphere_radius);
    double rho_t = ((param::rho_initial) / rho0) * rho_wd;

    double rho = b.getDensity();
    nominator += SQ((rho - rho_t) / rho_t);
  } // Loop over bodies

  mpi_utils::reduce_sum(nominator);
  mpi_utils::reduce_sum(nbs);

  double sdev = sqrt(nominator / (nbs - 1));
} // calculate standard deviation

/**
 * @brief      WVT displacement  dx
 *             [Diehl et al., PASA 2015]
 *
 *             dx = x + wvt_mu * wvt_acc
 *
 * @param      srch  The source's body holder
 */
void
wvt_displacement(body & source) {
  using namespace param;

  // Rename wvt_mu in case we want an iteration-dependent
  // value
  double wvt_mu_it = wvt_mu;
  double boundary = wvt_radius;

  if(physics::iteration > final_iteration) {
    int it = physics::iteration - final_iteration;
    double a = -0.1 * wvt_mu / (wvt_cool_down * 1.0);
    double b = 0.1 * wvt_mu;
    wvt_mu_it = std::max(a * it + b, 0.0);
  }

  // Check where particle will end up in next iteration
  point_t rp = source.coordinates() + wvt_mu_it * source.getAcceleration();

  double mass = source.mass();
  double r = magnitude(rp);

  // Always decrease "timestep" when particles move too far
  // out. This should be implemented globally eventually,
  // similar to a daptive timestep
  if(r / boundary > 1.01 && wvt_mu_it > 1e-10) {
    while(r / boundary > 1.01 && wvt_mu_it > 1e-10) {
      wvt_mu_it *= 0.5;
      rp = source.coordinates() + wvt_mu_it * source.getAcceleration();
      r = magnitude(rp);
    }
    // source.setAcceleration(0.0);
  }

  // Freeze particels in the outer edge
  if(boost::iequals(wvt_boundary, "frozen")) {
  }
  // Reflect particles on the sphere edge
  if(r / boundary >= 1.0) {
    point_t rs = cartesian_to_spherical(rp);
    rs[0] = rs[0] - 2.0 * (rs[0] - boundary);
    rp = spherical_to_cartesian(rs);
    // source.setAcceleration(0.0);
  }

  r = magnitude(rp);
  double rho0 = density_profiles::spherical_density_profile(0);
  double rho = (param::rho_initial) / rho0 *
               density_profiles::spherical_density_profile(r / sphere_radius);
  double h = sph_eta * kernels::kernel_width * pow(mass / rho, 1. / gdimension);

  // neighbor-based smoothing length calculation
  if(wvt_h_ngb) {
    if(gdimension == 1) {
      h = wvt_ngb * mass / (2.0 * rho);
    }
    else if(gdimension == 2) {
      h = sqrt(wvt_ngb * mass / (M_PI * rho));
    }
    else {
      h = pow(wvt_ngb * mass / (4.0 * M_PI * rho / 3.0), 1.0 / 3.0);
    }
  }
  if(rho == 0.0) {
    h = sphere_radius;
  }
  source.set_radius(h);
  source.set_coordinates(rp);
  source.setDensity(rho);
} // wvt_displacement

/**
 * @brief      Sets particle density to density profile
 *
 * @param      srch  The source's body holder
 */
void
wvt_set_density(body & source) {
  using namespace param;
  point_t rp = source.coordinates();
  double mass = source.mass();
  double r = magnitude(rp);

  double rho0 = density_profiles::spherical_density_profile(0);
  double rho = (param::rho_initial) / rho0 *
               density_profiles::spherical_density_profile(r / sphere_radius);
  double h = sph_eta * kernels::kernel_width * pow(mass / rho, 1. / gdimension);
  source.set_radius(h);
  source.setDensity(rho);
}

// wvt types and pointers
typedef void (*compute_quantity_t)(body &, std::vector<body *> &);
compute_quantity_t wvt_acceleration = wvt_acceleration_diehl;

/**
 * @brief      WVT method selector
 */
void
select() {
  using namespace param;
  if(boost::iequals(wvt_method, "diehl")) {
    wvt_acceleration = wvt_acceleration_diehl;
  }
  else if(boost::iequals(wvt_method, "arth")) {
    wvt_acceleration = wvt_acceleration_arth;
  }
  else {
    logm(error) << "ERROR: No WVT method specified";
    exit(2);
  }
} // select()

}; // namespace wvt

#endif // _wvt_h_
