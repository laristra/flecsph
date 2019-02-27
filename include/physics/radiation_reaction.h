/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
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
 * @file radiation_reaction.h
 * @brief Implementation of utilities for GW 
 *        radiation reaction
 */

#if 0
// TODO : Finish this Oct.6.2018

#ifndef _radiation_reaction_h_
#define _radiation_reaction_h_

#include <vector>

#include "params.h"
#include "utils.h"
#include "kernels.h"
#include "tree.h"

namespace radiation_reaction{

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

  double get_vec_dist(const point_t &p1, const point_t &p2) {
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

}; //radiation reaction

void binary_system_prop_pre_comp(){
}

void compute_gw_acc_part(){
}

void get_part_offset(){
}

int get_star_index(){
}

void get_ang_mom_vec(){
}

#endif // _radiation_reaction_h

#endif // Block comment
