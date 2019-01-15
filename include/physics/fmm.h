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
 * @file fmm.h
 * @brief Functions used in the FMM computation
 */

 #ifndef _fmm_h_
 #define _fmm_h_

 #include "params.h"
 #include "tree.h"

namespace fmm {
  using namespace param;

  /*
  * @brief Compute the gravitation interation
  * Return the computed value if needed for direct particle interaction
  *
  * The Sink is the one on which I compute the fc, dfcdr, dfcdrdr
  * The Source is the distant box
  */
  point_t gravitation_fc(
    point_t & fc,
    const point_t& sink_coordinates,
    const point_t& source_coordinates,
    const double& source_mass
  )
  {
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    point_t res = {};
    if(dist > 0){
      res = source_mass/(dist*dist*dist)*(source_coordinates-sink_coordinates);
    }
    fc = fc - res;
    return res;
  }

  /*
  * @brief Compute the Jacobian (dfcdr) matrix on the Sink from the Source
  */
  void gravitation_dfcdr(
    double dfcdr[9],
    const point_t& sink_coordinates,
    const point_t& source_coordinates,
    const double& source_mass
  )
  {
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    assert(dist > 0.0);
    point_t diffPos =  sink_coordinates - source_coordinates;
    double jacobicoeff = -source_mass/(dist*dist*dist);
    for(size_t i=0;i<gdimension;++i){
      for(size_t j=0;j<gdimension;++j){
        double valjacobian = 0.;
        if(i==j){
          valjacobian = jacobicoeff*(1-3*diffPos[i]*diffPos[j]/(dist*dist));
        }else{
          valjacobian = jacobicoeff*(-3*diffPos[i]*diffPos[j]/(dist*dist));
        }
        assert(!std::isnan(valjacobian));
        dfcdr[i*gdimension+j] += valjacobian;
      }
    }
  }

  /*
  * @brief Compute the Hessian (dfcdrdr) matrix on the Sink from the Source
  */
  void gravitation_dfcdrdr(
    double dfcdrdr[27],
    const point_t& sink_coordinates,
    const point_t& source_coordinates,
    const double& source_mass
  )
  {
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    assert(dist > 0.0);
    point_t diffPos =  sink_coordinates - source_coordinates;
    double hessiancoeff = -3.0*source_mass/(dist*dist*dist*dist*dist);
    for(size_t i=0;i<gdimension;++i){
      int matrixPos = i*gdimension*gdimension;
      for(size_t j=0;j<gdimension;++j){
        for(size_t k=0;k<gdimension;++k){
          int position = matrixPos+j*gdimension+k;
          double firstterm = 0.0;
          if(i==j){
            firstterm += diffPos[k];
          } // if
          if(j==k){
            firstterm += diffPos[i];
          } // if
          if(k==i){
            firstterm += diffPos[j];
          } // if
          double valhessian = hessiancoeff *
            ( 5.0/(dist*dist)*diffPos[i]*diffPos[j]*diffPos[k] - firstterm) ;
          dfcdrdr[position] += valhessian;
        } // for
      } // for
    } // for
  }

  /*
  * @brief Taylor expansion of degree 2 using the computed function, jacobi,
  * hessian and the targeted particle
  */
  void interation_c2p(
    double& fc,
    double dfcdr[9],
    double dfcdrdr[27],
    point_t cofm_coordinates,
    body_holder* sink
  )
  {
    point_t part_coordinates = sink->getBody()->coordinates();

    point_t diffPos = part_coordinates - cofm_coordinates;
    point_t grav = fc;
    // The Jacobi
    for(size_t i=0;i<gdimension;++i){
      for(size_t  j=0;j<gdimension;++j){
        grav[i] += dfcdr[i*gdimension+j]*diffPos[j];
      } // for
    } // for
    // The hessian
    double tmpMatrix[gdimension*gdimension] = {};
    for(size_t i=0;i<gdimension;++i){
      for(size_t j=0;j<gdimension;++j){
        for(size_t k=0;k<gdimension;++k){
          tmpMatrix[i*gdimension+j] +=
            diffPos[k]*dfcdrdr[i*gdimension*gdimension+j*gdimension+k];
        } // for
      } // for
    } // for
    double tmpVector[gdimension] = {};
    for(size_t i=0;i<gdimension;++i){
      for(size_t j=0;j<gdimension;++j){
        tmpVector[j] += tmpMatrix[i*gdimension+j]*diffPos[i];
      } // for
    } // for
    for(size_t i=0;i<gdimension;++i){
      grav[i] += 0.5*tmpVector[i];
    } // for
    sink->getBody()->setAcceleration(grav+sink->getBody()->getAcceleration());
  }

} // namespace fmm

#endif
