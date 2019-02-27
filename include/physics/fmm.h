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
  inline
  point_t gravitation_fc(
    point_t & fc,
    const point_t& sink_coordinates,
    const point_t& source_coordinates,
    const double& source_mass
  )
  {
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    point_t res = -source_mass/(dist*dist*dist)*
      (sink_coordinates-source_coordinates);
    fc += res;
    return res;
  }

  /*
  * @brief Compute the Jacobian (dfcdr) matrix on the Sink from the Source
  */
  inline
  void gravitation_dfcdr(
    double dfcdr[9],
    const point_t& sink_coordinates,
    const point_t& source_coordinates,
    const double& source_mass
  )
  {
    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    double dist_2 = dist*dist;
    point_t diffPos =  sink_coordinates - source_coordinates;
    double jacobicoeff = -source_mass/(dist_2*dist);
    for(int i = 0; i < 9; ++i){
      int a = i/3; int b = i%3;
      double valjacobian = jacobicoeff*((a==b)-3*diffPos[a]*diffPos[b]/(dist_2));
      dfcdr[i] += valjacobian;
    }
  }

  /*
  * @brief Compute the Hessian (dfcdrdr) matrix on the Sink from the Source
  */
  inline
  void gravitation_dfcdrdr(
    double dfcdrdr[27],
    const point_t& sink_coordinates,
    const point_t& source_coordinates,
    const double& source_mass
  )
  {

    double dist = flecsi::distance(sink_coordinates,source_coordinates);
    double dist_2 = dist*dist;
    point_t diffPos =  sink_coordinates - source_coordinates;
    double hessiancoeff = -3.0*source_mass/(dist_2*dist_2*dist);
    for(int i = 0 ; i < 27 ; ++i){
      int a = i/9; int b = (i%9)/3; int c = i%3;
      double term_1 = (a==b)*diffPos[c]+(c==a)*diffPos[b]+(b==c)*diffPos[a];
      double valhessian = hessiancoeff *
        ( 5.0/(dist_2)*diffPos[a]*diffPos[b]*diffPos[c] - term_1) ;
      dfcdrdr[i] += valhessian;
    }

/*
    for(size_t i=0;i<gdimension;++i){
      size_t matrixPos = i*gdimension*gdimension;
      for(size_t j=0;j<gdimension;++j){
        for(size_t k=0;k<gdimension;++k){
          size_t position = matrixPos+j*gdimension+k;
          double term_1 = (i==j)*diffPos[k]+(j==k)*diffPos[i]+(k==i)*diffPos[j];
          double valhessian = hessiancoeff *
            ( 5.0/(dist_2)*diffPos[i]*diffPos[j]*diffPos[k] - term_1) ;
          dfcdrdr[position] += valhessian;
        } // for
      } // for
    } // for*/
  }

  /*
  * @brief Taylor expansion of degree 2 using the computed function, jacobi,
  * hessian and the targeted particle
  */
  void interation_c2p(
    point_t& fc,
    double dfcdr[9],
    double dfcdrdr[27],
    point_t cofm_coordinates,
    body* sink
  )
  {
    point_t part_coordinates = sink->coordinates();

    point_t diffPos = part_coordinates - cofm_coordinates;
    point_t grav = fc;
    // The Jacobi
    for(size_t i=0;i<gdimension;++i){
      for(size_t j=0;j<gdimension;++j){
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
    sink->setAcceleration(grav+sink->getAcceleration());
  }

} // namespace fmm

#endif
