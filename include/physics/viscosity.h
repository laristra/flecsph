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

#ifndef _viscosity_h_
#define _viscosity_h_

#include <vector>

#include <boost/algorithm/string.hpp>

#define SQ(x) ((x)*(x))
#define QU(x) ((x)*(x)*(x)*(x))

namespace viscosity{
  using namespace param;
  static const double TINY = 1e10*DBL_MIN;



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
  mu(
      const double & h_ab,
      const space_vector_t & vel_ab,
      const space_vector_t & pos_ab)
  {

    using namespace param;
    double result = 0.0;
    double dotproduct = flecsi::dot(vel_ab, pos_ab);
    double dist2 = flecsi::dot(pos_ab,pos_ab);
    result = h_ab*dotproduct / (dist2 + sph_viscosity_epsilon*h_ab*h_ab + TINY);

    //mpi_assert(result < 0.0);
    return result*(dotproduct < 0.0);
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
  artificial_viscosity(
    const double & alpha_a,
    const double & alpha_b,
    const double & rho_ab,
    const double & c_ab,
    const double & mu_ab)
  {
    using namespace param;
    double res = ( -sph_viscosity_alpha*c_ab
                  + sph_viscosity_beta*mu_ab)*mu_ab/rho_ab;
    //mpi_assert(res>=0.0);
    return res;
  }

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
  viscosity_cullen(
    const double & alpha_a,
    const double & alpha_b,
    const double & rho_ab,
    const double & c_ab,
    const double & mu_ab)
  {
    using namespace param;
    double alpha = (alpha_a+alpha_b)/2.0;
    double res = ( -alpha*c_ab
                  + 2.0*alpha*mu_ab)*mu_ab/rho_ab;
    //mpi_assert(res>=0.0);
    return res;
  }

  typedef double (*viscosity_function_t)(const body &, const body &);
  viscosity_function_t viscosity = nullptr;

  /**
   * @brief Viscosity selector
   * @param kstr Viscosity string descriptor
   */
  void select(const std::string& kstr)
  {
    if (boost::iequals(kstr,"artificial_viscosity")){
      viscosity = nullptr;
    } else if (boost::iequals(kstr,"artificial_cullen")){
      viscosity = nullptr;
    } else{
      clog_fatal("Bad viscosity parameter"<<std::endl);
    }
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
   * @uses       alpha_max  global parameter
   */
  void
  compute_alpha(
    body& particle,
    std::vector<body*>& nbs)
  {
    using namespace param;
    using namespace kernels;
    // this particle (index 'a')
    const double c_a = particle.getSoundspeed(),
                 h_a = particle.radius(),
             alpha_a = particle.getAlpha();
    const point_t pos_a = particle.coordinates(),
                    v_a = particle.getVelocityhalf();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double c_a[n_nb];
    point_t pos_a[n_nb], v_a[n_nb], ;

    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      pos_a[b] = (pos_a - nb->coordinates())/flecsi::distance(pos_a, pos_b);
      v_a[b]   = v_a - nb->getVelocity();
      c_a[b] = 0.5*(c_a + nb->getSoundspeed());
    }

    // compute signal velocity
    double vsig = 0.0;
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double temp = std::min(0.0,flecsi::dot(v_a[b], pos_a[b]));
      temp = c_a[b] - temp;
      if (temp > vsig){
        vsig = temp;
      }
    }
    double Atrig = A_trigger(particle, nbs);
    double alpha_loc = alpha_max*Atrig/(Atrig + SQ(vsig)/SQ(h_a));

    if (alpha_a < alpha_loc){
      particle.setAlpha(alpha_loc);
    }
    else if (alpha_a > alpha_loc){
      double decayt = h_a/(2.0*viscosity_l*vsig);
      double dalphadt = (alpha_loc - alpha_a)/decayt;
      particle.setAlpha(alpha_a + dalphadt*physics::dt);
    }
    // compute the final answer
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double h_ab = .5*(h_a + h_[b]);
      DiWab = sph_kernel_gradient(pos_a - pos_[b],h_ab);
      result += m_[b]*flecsi::dot(v_a[b],DiWab);
    }
    result /= rho_a;
    particle.setDivergenceV(result);
  } // compute_alpha

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
  A_trigger(
    body& particle,
    std::vector<body*>& nbs)
  {

    using namespace param;
    const double DivV_a_old = particle.getDivergenceV();
    double dDivVdt = 0.0;
    double result = 0.0;
    double xi = 0.0;

    compute_DivergenceV(particle, nbs);
    dDivVdt = (particle.getDivergenceV() - DivV_old)/physics::dt;

    result = std::max(-dDivVdt,0.0);

    xi = compute_xi(particle,nbs);

    result = xi*result;
    return result;
  } // A_trigger

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
                  rho_a = particle.getDensity();
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
      m_[b]    = nb->mass() * (pos_[b]!=pos_a); // if same particle, m_b->0
    }

    // precompute velocity difference and kernel gradients
    // compute R_a
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      const space_vector_t pos_ab = point_to_vector(pos_a - pos_[b]);
      double h_ab = .5*(h_a + h_[b]);
      v_a_[b]  = point_to_vector(v_[b] - v_a);
      DiWa_[b] = sph_kernel_gradient(pos_ab,h_ab);

      double Wab =  sph_kernel_function(flecsi::distance(pos_a, pos_[b]),h_ab);
      R_a += signnum_c(divV_[b])*m_[b]*Wab;
    }

    // calculate the gradient of velocity matrix
    for(int i = 0; i < gdimension; i++){
      for(int j = 0; j < gdimension; j++){
        for(int b = 0 ; b < n_nb; ++b){
          gradV_a[(gdimension*i)+j] += m_[b]*v_a_[b][i]*DiWa_[b][j];
        }
      }
    }
    gradV_a /= rho_a;

    // traceless symmetric part of velocity gradient
    for(int i = 0; i < gdimension; i++){
      for(int j = 0; j < gdimension; j++){
        SymT_a[(gdimension*i)+j] = gradV_a[(gdimension*i)+j]+gradV_a[(gdimension*j)+i];
      }
    }
    SymT_a /= 2.0;
    for(int i = 0; i < gdimension; i++){
      SymT_a[(gdimension*i)+i] -= divV_a/gdimension;
    }

    // trace of (S S^dagger)
    for(int i = 0; i < gdimension; i++){
      for(int j = 0; j < gdimension; j++){
        traceSS_a += SymT_a[(gdimension*i)+j]*SymT_a[(gdimension*i)+j];
      }
    }

    // compute the final answer
    result = SQ(2.0*QU(1.0-R_a)*divV_a);
    return result = result/(result + traceSS_a);
  } // compute_xi

  /**
   * @brief      Computes the divergence of velocity for a particle a
   *
   * @param      particle  The particle body
   * @param      nbs       Vector of neighbor particles
   */
  void
  compute_DivergenceV(
      body& particle,
      std::vector<body*>& nbs)
  {
    using namespace param;
    using namespace kernels;
    // this particle (index 'a')
    const double h_a = particle.radius(),
               rho_a = particle.getDensity();
    const point_t pos_a = particle.coordinates(),
                  v_a = particle.getVelocityhalf();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double h_[n_nb],m_[n_nb];
    point_t pos_[n_nb], v_a[n_nb], DiWab;

    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      pos_[b] = nb->coordinates();
      v_a[b]  = nb->getVelocity() - v_a;
      h_[b]   = nb->radius();
      m_[b]   = nb->mass() * (pos_[b]!=pos_a); // if same particle, m_b->0
    }

    // compute the final answer
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double h_ab = .5*(h_a + h_[b]);
      DiWab = sph_kernel_gradient(pos_a - pos_[b],h_ab);
      result += m_[b]*flecsi::dot(v_a[b],DiWab);
    }
    result /= rho_a;
    particle.setDivergenceV(result);

  } // compute_DivergenceV

  /**
   * @brief      return the sign of double given to function
   *
   * @param      x          The double
   *
   * @return     sign of double
   *
   */
  double signnum_c(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return x;
  }
}; // viscosity
#undef SQ
#undef QU

#endif // _viscosity_h_
