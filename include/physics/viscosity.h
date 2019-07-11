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

  typedef double (*viscosity_function_t)(const body &, const body &);
  viscosity_function_t viscosity = nullptr; //artificial_viscosity;

  /**
   * @brief Viscosity selector
   * @param kstr Viscosity string descriptor
   */
  void select(const std::string& kstr)
  {
    if (boost::iequals(kstr,"artificial_viscosity")){
      viscosity = nullptr; //artificial_viscosity;
    }else{
      clog_fatal("Bad viscosity parameter"<<std::endl);
    }
  }

  /**
   * @brief      alpha_loc_i for the artificial viscosity
   * From Cullen'10 (arXiv:1006.1524) -
   * Inviscid SPH, eq.(14)
   *
   * @param      srch       The source particle
   * @param      nbsh       The neighbor particle
   *
   * @return     Contribution for alpha_loc_i
   *
   * @uses       alpha_max  global parameter
   */
  inline double
  alpha_loc(
      )
  {

    using namespace param;
    double result = 0.0;

    double dotproduct = alpha_max*flecsi::dot(vel_ab, pos_ab);
    double dist2 = flecsi::dot(pos_ab,pos_ab);
    result = h_ab*dotproduct / (dist2 + sph_viscosity_epsilon*h_ab*h_ab + TINY);

    //mpi_assert(result < 0.0);
    return result*(dotproduct < 0.0);
  } // alpha_loc

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
    const double DivV_old = particle.getDivergenceV();
    double dDivVdt = 0.0;
    double result = 0.0;
    double eeeeeee = 0.0;

    compute_DivergenceV(particle);
    dDivVdt = (particle.getDivergenceV()-DivV_old)/physics::dt;

    result = std::max(dDivVdt,0.0);

    eeeeeee = 

    double dotproduct = alpha_max*flecsi::dot(vel_ab, pos_ab);
    double dist2 = flecsi::dot(pos_ab,pos_ab);
    result = h_ab*dotproduct / (dist2 + sph_viscosity_epsilon*h_ab*h_ab + TINY);

    //mpi_assert(result < 0.0);
    return result*(dotproduct < 0.0);
  } // alpha_loc

}; // viscosity

#endif // _viscosity_h_
