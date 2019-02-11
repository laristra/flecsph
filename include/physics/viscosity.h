/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
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
      const double & h_a,
      const double & h_b,
      const point_t & vel_a,
      const point_t & vel_b,
      const point_t & pos_a,
      const point_t & pos_b)
  {

    using namespace param;
    double result = 0.0;
    double h_ab = .5*(h_a + h_b);
    space_vector_t vecVelocity = flecsi::point_to_vector(vel_a - vel_b);
    space_vector_t vecPosition = flecsi::point_to_vector(pos_a - pos_b);
    double dotproduct = flecsi::dot(vecVelocity,vecPosition);

    double dist = flecsi::distance(pos_a, pos_b);
    result = h_ab*dotproduct / (dist*dist + sph_viscosity_epsilon*h_ab*h_ab);

    //mpi_assert(result < 0.0);
    return result*(dotproduct >= 0.0);
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
    const double & rho_a, 
    const double & rho_b,
    const double & c_a,
    const double & c_b,
    const double & mu_ab)
  {
    using namespace param;
    double rho_ab = .5*(rho_a + rho_b);
    double c_ab = .5*(c_a + c_b);
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

}; // viscosity

#endif // _viscosity_h_
