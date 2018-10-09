/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@ @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@///// /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@      /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@ /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////  /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@      /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@      /@@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file analysis.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic analysis
 */

#ifndef _PHYSICS_ANALYSIS_H_
#define _PHYSICS_ANALYSIS_H_

#include <vector>
#include "params.h"

//#include "physics.h"

namespace analysis{

  point_t linear_momentum;
  point_t linear_velocity;
  point_t total_ang_mom;
  point_t part_mom;
  point_t part_position;
  double total_mass;
  double total_energy;
  double velocity_part;

  /**
   * @brief      Compute the linear momentum
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_lin_momentum(
      std::vector<body_holder*>& bodies)
  {
    linear_momentum = {0};
    for(auto nbh: bodies) {
      linear_momentum += nbh->getBody()->getLinMomentum();
    }
    reduce_sum(linear_momentum);
  }

  /**
   * @brief      Compute the total mass
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_mass(
      std::vector<body_holder*>& bodies)
  {
    total_mass = 0.;
    for(auto nbh: bodies) {
      total_mass += nbh->getBody()->getMass();
    }
    reduce_sum(total_mass);
  }

  /**
   * @brief      Compute the total energy (internal + kinetic)
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_energy(
      std::vector<body_holder*>& bodies)
  {
    using namespace param;
    total_energy = 0.;
    if (thermokinetic_formulation) {
      for(auto nbh: bodies)
        total_energy += nbh->getBody()->getMass()*nbh->getBody()->getTotalenergy();
    }
    else {
      for(auto nbh: bodies) {
        total_energy += nbh->getBody()->getMass()*nbh->getBody()->getInternalenergy();
        linear_velocity = nbh->getBody()->getVelocity();
        velocity_part = 0.;
        for(size_t i = 0 ; i < gdimension ; ++i){
          velocity_part += pow(linear_velocity[i],2);
          part_position = nbh->getBody()->getPosition();
        }
        total_energy += 1./2.*velocity_part*nbh->getBody()->getMass();
      }
    }
    reduce_sum(total_energy);
  }

  /**
   * @brief      Compute the total energy (internal + kinetic)
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_ang_mom(
      std::vector<body_holder*>& bodies)
  {
    double part_ang_x;
    double part_ang_y;
    double part_ang_z;
    total_ang_mom = {0};
    part_mom = {0};
    part_position = {0};
    for(auto nbh: bodies) {
      part_mom = nbh->getBody()->getLinMomentum();
      part_position = nbh->getBody()->getPosition();
      if(gdimension==3){
        part_ang_x = part_position[1]*part_mom[2]-part_position[2]*part_mom[1];
        part_ang_y = -part_position[0]*part_mom[2]+part_position[2]*part_mom[0];
        part_ang_z = part_position[0]*part_mom[1]-part_position[1]*part_mom[0];
        total_ang_mom += {part_ang_x,part_ang_y,part_ang_z};
      } else if(gdimension==2){
        part_ang_z = part_position[0]*part_mom[1]-part_position[1]*part_mom[0];
        total_ang_mom += part_ang_z;
      } else if(gdimension==1){
        total_ang_mom += 0.;
      }
    }
    reduce_sum(total_ang_mom);
  }

  /**
   * @brief Rolling screen output
   */
  void
  screen_output(int rank)
  {
    using namespace param;
    static int count = 0;
    const int screen_length = 40;
    if (out_screen_every > 0 || physics::iteration % out_screen_every == 0) {
      (++count-1)%screen_length  ||
      rank || clog(info)<< "#-- iteration:               time:" <<std::endl;
      rank || clog(info)
        << std::setw(14) << physics::iteration
        << std::setw(20) << std::scientific << std::setprecision(12)
        << physics::totaltime << std::endl;
    }
  }


  /**
   * @brief Periodic file output
   *
   * Outputs all scalar reductions as single formatted line to a file. When
   * creating the file, writes a header to indicate which quantities are
   * output in which column (because their order and quantity may change
   * between revisions)
   * E.g.:
   * -- >> example output file >> -----------------------------------------
   * # Scalar reductions:
   * # 1:iteration 2:time 3:energy 4:mom_x 5:mom_y 6:mom_z
   * 0  0.0   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * 10 0.1   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * 20 0.2   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * 30 0.3   1.0000000e+03  0.00000000e+00  0.00000000e+00  0.00000000e+00
   * ...
   * -- << end output file <<<< -------------------------------------------
   */
  void
  scalar_output(const char * filename = NULL)
  {
    static bool first_time = true;

    // output only from rank #0
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (rank != 0) return;

    if (first_time) {
      // Generate and output the header
      std::ostringstream oss_header;
      oss_header
        << "# Scalar reductions: " <<std::endl
        << "# 1:iteration 2:time 3:mass 4:total_energy 5:mom_x ";

      // The momentum depends on dimension
      if(gdimension > 1){
        oss_header<<"6:mom_y ";
      }
      if(gdimension == 3){
        oss_header<<"7:mom_z ";
      }
      // The angular momentum depends on dimension
      if(gdimension == 2){
        oss_header<<"7:ang_mom_z";
      }
      if(gdimension == 3){
        oss_header<<"8:ang_mom_x 9:ang_mom_y 10:ang_mom_y ";
      }
      oss_header<<std::endl;

      std::ofstream out(filename);
      out << oss_header.str();
      out.close();

      first_time = false;
    }

    std::ostringstream oss_data;
    oss_data << std::setw(14) << physics::iteration
      << std::setw(20) << std::scientific << std::setprecision(12)
      << physics::totaltime << std::setw(20) << total_mass
      << std::setw(20) << total_energy;
    for(size_t i = 0 ; i < gdimension ; ++i){
      oss_data
        << std::setw(20) << std::scientific << std::setprecision(12)
        << linear_momentum[i] <<" ";
    }
    if(gdimension==2){
      oss_data
        << std::setw(20) << std::scientific << std::setprecision(12)
        << total_ang_mom[0] <<" ";
    }
    if(gdimension==3){
      for(size_t i = 0 ; i < gdimension ; ++i){
        oss_data
          << std::setw(20) << std::scientific << std::setprecision(12)
          << total_ang_mom[i] <<" ";
      }
    }
    oss_data << std::endl;

    // Open file in append mode
    std::ofstream out(filename,std::ios_base::app);
    out << oss_data.str();
    out.close();

  } // scalar output

}; // physics

#endif // _PHYSICS_ANALYSIS_H_
