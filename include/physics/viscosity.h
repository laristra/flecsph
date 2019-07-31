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
    double alpha = 0.5*(alpha_a+alpha_b);
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
   * @uses       sph_viscosity_alpha_max  global parameter
   */
  void
  initialize_alpha(
    body& particle)
  {
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
    particle.setTraceSS(traceSS_a);
    // compute the final answer
    result = SQ(2.0*QU(1.0-R_a)*divV_a);
    if(result == 0)
      return 0;
    return result = result/(result + traceSS_a);
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
  A_trigger(
    body& particle,
    std::vector<body*>& nbs,
    const double& DivV_a_new)
  {

    using namespace param;
    const double DivV_a_old = particle.getDivergenceV();
    double dDivVdt = 0.0;
    double result = 0.0;
    double xi = 0.0;

    //compute_DivergenceV(particle, nbs);
    dDivVdt = (DivV_a_new - DivV_a_old)/physics::dt;

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
  compute_alpha(
    body& particle,
    std::vector<body*>& nbs)
  {
    using namespace param;
    using namespace kernels;
    // this particle (index 'a')
    const double c_a = particle.getSoundspeed(),
                 h_a = particle.radius(),
             alpha_a = particle.getAlpha(),
               rho_a = particle.getDensity();
    const point_t pos_a = particle.coordinates(),
                    v_a = particle.getVelocity();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double h_[n_nb], m_[n_nb], c_a_[n_nb];
    point_t pos_[n_nb], pos_a_[n_nb], v_a_[n_nb], DiWab;

    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      pos_[b]   = nb->coordinates();
      pos_a_[b] = (pos_a - pos_[b])/flecsi::distance(pos_a, pos_[b]);
      v_a_[b]   = v_a - nb->getVelocity();
      c_a_[b] = 0.5*(c_a + nb->getSoundspeed());
      h_[b]     = nb->radius();
      m_[b]     = nb->mass();
    }

    // compute signal velocity
    double vsig = 0.0;
    for(int b = 0 ; b < n_nb; ++b){
      double dotval = v_a_[b][0]*pos_a_[b][0];
      for (unsigned short i=1; i<gdimension; ++i)
        dotval += v_a_[b][i]*pos_a_[b][i];
      double temp = std::min(0.0,dotval);
      temp = c_a_[b] - temp;

      if (temp > vsig){
        vsig = temp;
      }
    }
    double result = 0.0;
    // compute the divergence
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double h_ab = .5*(h_a + h_[b]);
      DiWab = sph_kernel_gradient(pos_a - pos_[b],h_ab);
      double dotval = -v_a_[b][0]*DiWab[0];
      for (unsigned short i=1; i<gdimension; ++i)
        dotval += -v_a_[b][i]*DiWab[i];
      result += m_[b]*dotval;
    }
    result /= rho_a;
    double Atrig = A_trigger(particle, nbs, result);
    double alpha_loc = sph_viscosity_alpha_max*Atrig/(Atrig + SQ(vsig)/SQ(h_a));

    if (alpha_a <= alpha_loc){
      particle.setAlpha(alpha_loc);
    }
    else {
      double decayt = h_a/(2.0*sph_viscosity_l*vsig);
      double dalphadt = (alpha_loc - alpha_a)/decayt;
      particle.setAlpha(alpha_a + dalphadt*physics::dt);
    }

    particle.setDivergenceV(result);
  } // compute_alpha


}; // viscosity
#undef SQ
#undef QU

#endif // _viscosity_h_
