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
 * @file viscosity.h
 * @author Julien Loiseau
 * @date October 2018
 * @brief Viscosities implementations
 */

#pragma once

#include <vector>
#include <boost/algorithm/string.hpp>

#define SQ(x) ((x)*(x))
#define QU(x) ((x)*(x)*(x)*(x))

namespace viscosity {
using namespace param;
static const double TINY = 1e10*DBL_MIN;

// Generic template: artificial viscosity function
template<param::sph_viscosity_keyword K>
double
viscosity_function(
  const double alpha_ab,
  const double rho_ab,
  const double c_ab,
  const double mu_ab);

// Function pointers
typedef double (*viscosity_function_t)(
  const double alpha_ab,
  const double rho_ab,
  const double c_ab,
  const double mu_ab);

/**
 * @brief      return the sign of double given to function
 *
 * @param      x          The double
 *
 * @return     sign of double
 *
 */
template<typename T>
T signnum_c(T x) {
  if (x > 0.0) return T(1);
  if (x < 0.0) return T(-1);
  return x;
}

/**
 * @brief      mu_ij for the artificial viscosity
 * From Rosswog'09 (arXiv:0903.5075) -
 * Astrophysical Smoothed Particle Hydrodynamics, eq.(60)
 *
 * @param      srch     The source particle
 * @param      nbsh     The neighbor particle
 *
 * @return     Contribution for mu_ij of this neighbor
 *
 * @uses       epsilon  global parameter
 */
inline double
mu(const double & h_ab, const point_t & vel_ab, const point_t & pos_ab) {

  using namespace param;
  double result = 0.0;
  double dotproduct = flecsi::dot(vel_ab, pos_ab);
  double dist2 = flecsi::dot(pos_ab, pos_ab);
  result =
    h_ab * dotproduct / (dist2 + sph_viscosity_epsilon * h_ab * h_ab + TINY);

  // mpi_assert(result < 0.0);
  return result * (dotproduct < 0.0);
} // mu

/**
 * @brief      Artificial viscosity term, Pi_ab
 * From Rosswog'09 (arXiv:0903.5075) -
 * Astrophysical Smoothed Particle Hydrodynamics, eq.(59)
 *
 * @param      srch  The source particle
 * @param      nbsh  The neighbor particle
 *
 * @return     The artificial viscosity contribution
 */
inline double
artificial_viscosity(const double & rho_ab,
  const double & c_ab,
  const double & mu_ab) {
  using namespace param;
  double res =
    (-sph_viscosity_alpha * c_ab + sph_viscosity_beta * mu_ab) * mu_ab / rho_ab;
  // mpi_assert(res>=0.0);
  return res;
}

/**
 * @brief      Artificial viscosity term, Pi_ab
 * From Rosswog'09 (arXiv:0903.5075) -
 * Astrophysical Smoothed Particle Hydrodynamics, eq.(59)
 *
 * @param      alpha_ab  not used
 * @param      rho_ab    average density
 * @param      c_ab      average speed of sound
 * @param      mu_ab     average mu (see the mu function)
 *
 * @return     The artificial viscosity contribution
 */
template<> double
viscosity_function<param::visc_constant> (
  const double alpha_ab,  // this parameter is ignored
  const double rho_ab,
  const double c_ab,
  const double mu_ab)
{
  using namespace param;
  return ( -sph_viscosity_alpha*c_ab + sph_viscosity_beta*mu_ab)*mu_ab/rho_ab;
}

/**
 * @brief      Artificial viscosity term, Pi_ab
 * From Cullen & Dehnen (2010) (arXiv:1006.1524), "Inviscid SPH", eq.(4)
 *
 * @param      alpha_ab  average viscosity-alpha parameter
 * @param      rho_ab    average density
 * @param      c_ab      average speed of sound
 * @param      mu_ab     average mu (see the mu function)
 *
 * @return     The artificial viscosity contribution
 */
template<> double
viscosity_function<param::visc_cullen> (const double alpha_ab,
  const double rho_ab,
  const double c_ab,
  const double mu_ab) {
  using namespace param;
  return -alpha_ab*(c_ab - 2.0*mu_ab)*mu_ab/rho_ab;
}

/**
 * @brief      alpha parameter for the artificial viscosity
 * From Cullen'10 (arXiv:1006.1524) -
 *
 * @param      srch       The source particle
 * @param      nbsh       The neighbor particle
 *
 * @return
 *
 * @uses       sph_viscosity        global parameter
 * @uses       sph_viscosity_alpha  global parameter
 */
void
initialize_alpha(
  body& particle)
{
  if (param::sph_viscosity == param::visc_constant)
    particle.setAlpha(sph_viscosity_alpha);
  else
    particle.setAlpha(0.0);
} // initialize_alpha

/**
 * @brief      xi_a for the artificial viscosity:
 *             calc R_a
 * From Cullen'10 (arXiv:1006.1524) -
 * Inviscid SPH, eqs.(17, 18)
 *
 * @param      srch       The source particle
 * @param      nbsh       The neighbor particle
 *
 * @return     limiter to reduce unwanted dissipation
 *
 */
inline double
compute_xi(
  body& particle,
  std::vector<body*>& nbs)
{
  using namespace param;
  using namespace kernels;

  // this particle (index 'a')
  const double divV_a = particle.getDivergenceV(),
                rho_a = particle.getDensity(),
                  h_a = particle.radius();
  const point_t pos_a = particle.coordinates(),
                  v_a = particle.getVelocity();
  double gradV_a[gdimension*gdimension];
  double SymT_a[gdimension*gdimension];
  double traceSS_a = 0.0;
  double R_a = 0.0;
  double result = 0.0;

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  double h_[n_nb], m_[n_nb], divV_[n_nb];
  point_t pos_[n_nb], v_[n_nb], v_a_[n_nb], DiWa_[n_nb];

  for(int b = 0; b < n_nb; ++b) {
    const body * const nb = nbs[b];
    pos_[b]  = nb->coordinates();
    v_[b]    = nb->getVelocity();
    h_[b]    = nb->radius();
    divV_[b] = nb->getDivergenceV();
    m_[b]    = nb->mass();
  }

  // precompute velocity difference and kernel gradients
  // compute R_a
  for(int b = 0 ; b < n_nb; ++b){
    point_t pos_ab = pos_a - pos_[b];
    double h_ab = .5*(h_a + h_[b]);
    v_a_[b]  = v_a - v_[b];
    DiWa_[b] = sph_kernel_gradient(pos_ab,h_ab);

    double Wab =  sph_kernel_function(flecsi::distance(pos_a, pos_[b]),h_ab);
    R_a += signnum_c(divV_[b])*m_[b]*Wab;
  }
  R_a /= rho_a;

  // calculate the gradient of velocity matrix
  for(int i = 0; i < gdimension; i++){
    for(int j = 0; j < gdimension; j++){
      gradV_a[(gdimension*i)+j] = 0;
      for(int b = 0 ; b < n_nb; ++b){
        gradV_a[(gdimension*i)+j] += m_[b]*(-v_a_[b][i])*DiWa_[b][j];
      }
      gradV_a[(gdimension*i)+j] /= rho_a;
    }
  }
  particle.setGradV(gradV_a[0]);

  // traceless symmetric part of velocity gradient
  for(int i = 0; i < gdimension; i++){
    for(int j = 0; j < gdimension; j++){
      SymT_a[(gdimension*i)+j] = 0.5*(gradV_a[(gdimension*i)+j]+gradV_a[(gdimension*j)+i]);
    }
    SymT_a[(gdimension*i)+i] -= divV_a/gdimension;
    for(int j = 0; j < gdimension; j++){
      traceSS_a += SymT_a[(gdimension*i)+j]*SymT_a[(gdimension*i)+j];
    }
  }
  //traceSS_a = 0.0; // DEBUG
  particle.setTraceSS(traceSS_a);
  // compute the final answer
  result = SQ(2.0*QU(1.0-R_a)*divV_a);
  return (result < 1e-16) ? (0.0) : result/(result + traceSS_a);
} // compute_xi

/**
 * @brief      A_i for the artificial viscosity
 * From Cullen'10 (arXiv:1006.1524) -
 * Inviscid SPH, eq.(13)
 *
 * @param      srch       The source particle
 * @param      nbsh       The neighbor particle
 *
 * @return     Trigger for  viscosity
 *
 */
inline double
A_trigger( body& particle,
  std::vector<body*>& nbs,
  const double& DivV_a_new) {

  using namespace param;
  const double DivV_a_old = particle.getDivergenceV();
  double dDivVdt = 0.0;
  double result = 0.0;
  double xi = 0.0;

  dDivVdt = particle.getDdivvdt();

  ///result = std::max(-DivV_a_old,0.0); // DEBUG
  result = std::max(-dDivVdt,0.0);

  xi = compute_xi(particle,nbs);
  particle.setXi(xi);
  particle.setTrigger(result);
  result = xi*result;
  return result;
} // A_trigger

/**
 * @brief      alpha parameter for the artificial viscosity
 * From Cullen'10 (arXiv:1006.1524) -
 *
 * @param      srch       The source particle
 * @param      nbsh       The neighbor particle
 *
 * @return
 *
 * @uses       sph_viscosity_alpha_max  global parameter
 */
void
compute_alpha(body & particle, std::vector<body *> & nbs) {
  using namespace param;
  using namespace kernels;
  using namespace flecsi;
  // this particle (index 'a')
  const double c_a = particle.getSoundspeed(),
               h_a = particle.radius(),
           alpha_a = particle.getAlpha(),
             rho_a = particle.getDensity();
  const point_t pos_a = particle.coordinates(),
                  v_a = particle.getVelocity();

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  double h_[n_nb], m_[n_nb], c_a_[n_nb], rho_[n_nb];
  point_t pos_[n_nb], n_a_[n_nb], v_a_[n_nb], DiWab;

  for(int b = 0; b < n_nb; ++b) {
    const body * const nb = nbs[b];
    pos_[b]  = nb->coordinates();
    rho_[b]  = nb->getDensity();
    n_a_[b] = (pos_a - pos_[b])/distance(pos_a, pos_[b]);
    v_a_[b]   = v_a - nb->getVelocity();
    c_a_[b]   = std::max(c_a, nb->getSoundspeed());
    h_[b]    = nb->radius();
    m_[b]    = nb->mass();
  }

  // compute signal velocity
  double vsig = 0.0;
  for(int b = 0 ; b < n_nb; ++b){
    vsig = std::max(vsig, c_a_[b] - std::min(dot(v_a_[b],n_a_[b]),0.0));
  }

  double div_v = particle.getDivergenceV();
  double Atrig = A_trigger(particle, nbs, div_v);
  double alpha_loc = sph_viscosity_alpha_max 
                   * Atrig / (sph_viscosity_delta*SQ(vsig/h_a) + Atrig);

  if (alpha_a <= alpha_loc) {
    particle.setAlpha(alpha_loc);
  }
  else {
    double decayt = h_a/(2.0*sph_viscosity_l*vsig);
    double dalphadt = (alpha_a - alpha_loc)/decayt;
    particle.setAlpha(alpha_a*exp(-dalphadt*physics::dt));
  }

} // compute_alpha


#ifdef sph_viscosity
# define   sph_artificial_viscosity   viscosity_function<param::sph_viscosity>
#else
  viscosity_function_t sph_artificial_viscosity = nullptr;
#endif

/**
 * @brief Viscosity selector
 */
void select() {
#ifndef sph_viscosity
  using namespace param;
  switch(sph_viscosity) {
  case (visc_constant):
    sph_artificial_viscosity = viscosity_function<visc_constant>;
    break;
  case (visc_cullen):
    sph_artificial_viscosity = viscosity_function<visc_cullen>;
    break;
  default:
    log_fatal("Bad viscosity parameter" << std::endl);
  }
#endif
}


}; // namespace viscosity
#undef SQ
#undef QU

