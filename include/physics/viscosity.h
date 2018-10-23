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
  double 
  mu(
      body* source, 
      body* nb)
  {  
    using namespace param;
    double result = 0.0;
    double h_ij = .5*(source->getSmoothinglength()+nb->getSmoothinglength()); 
    space_vector_t vecVelocity = flecsi::point_to_vector(
        source->getVelocityhalf() - nb->getVelocityhalf());
    space_vector_t vecPosition = flecsi::point_to_vector(
        source->getPosition() - nb->getPosition());
    double dotproduct = flecsi::dot(vecVelocity,vecPosition);

    if(dotproduct >= 0.0)
      return result;
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    result = h_ij*dotproduct / (dist*dist + sph_viscosity_epsilon*h_ij*h_ij);
    
    mpi_assert(result < 0.0);
    return result; 
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
  double 
  artificial_viscosity(
    body* source, 
    body* nb)
  {
    using namespace param;
    double rho_ij = (1./2.)*(source->getDensity()+nb->getDensity());
    double c_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
    double mu_ij = mu(source,nb);
    double res = ( -sph_viscosity_alpha*c_ij*mu_ij
                  + sph_viscosity_beta*mu_ij*mu_ij)/rho_ij;
    mpi_assert(res>=0.0);
    return res;
  }

}; // viscosity

#endif // _viscosity_h_
