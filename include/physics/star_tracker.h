/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 * ~--------------------------------------------------------------------------~*/

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

/*
 * @file star_tracker.h
 * @brief Implementation of star tracker
 */

#ifndef _star_tracker_h_
#define _star_tracker_h_

#include <vector>

#include "params.h"
#include "utils.h"
#include "kernels.h"
#include "tree.h"

namespace star_tracker{

  //Need to: 
  //     find particle's max density
  //     find star's center of mass
  //     find whole system's center of mass
  //     Get angular momentum of particles

  //Struct for tracking positions of star system

  typedef struct {
      double mass; //Mass of star
      double radius; //Radius of star
      double com[gdimension]; //Center of mass of star
      double spin[gdimension]; //Angular momentue of star
  } singleStarData_t;


  //Struct for binary system
  #define NSTARS 2 //Define number of stars. Two is a default value

  typedef struct {
      double total_mass; //Total mass of system
      double reduced_mass; //Dimensionaless reduced mass
      double acc_com; //Acceleration per star
      double com[gdimension]; //Center of mass of star
      double spin[gdimension]; //Angular momentue of star
      double offset[NSTARS][gdimension]; //Offset; r_star - r_com
      double norm_offset[gdimension]; //Norm of offset
      double separation;
  } binaryStarData_t;

  double vec_dist(const point_t &p1, const point_t &p2) {
    if (gdimension==3) {
      return sqrt((p2[0]-p1[0])*(p2[0]-p1[0])
                  +(p2[1]-p1[1])*(p2[1]-p1[1])
                  +(p2[2]-p1[2])*(p2[2]-p1[2]));
    } else if (gdimension==2) {
      return sqrt((p2[0]-p1[0])*(p2[0]-p1[0])
                  +(p2[1]-p1[1])*(p2[1]-p1[1]));
    } else
      return abs(p2[0]-p1[0]);
  }

  void find_maxRho_parts(body_holder* srch) {
  }

  void find_starRho(body_holder* srch) {
  }

  void find_systemRho(body_holder* srch) {
  }

  void find_angMom_parts(body_holder* srch) {
  }

}; //star_tracker

#endif // _star_tracker_h
