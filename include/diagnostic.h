/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file diagnostic.h
 * @author Julien Loiseau
 * @date October 2018
 * @brief Generic particle diagnostic functions
 */

#ifndef _PHYSICS_DIAGNOSTIC_H_
#define _PHYSICS_DIAGNOSTIC_H_

#include <vector>
#include "params.h"

namespace diagnostic {

  uint64_t N_min, N_max, N_average;
  uint64_t N_ghosts;
  double h_min, h_max, h_average;
  double min_dist, average_dist_in_h;
  size_t id_N_min, id_N_max, id_h_min, id_h_max;
  double V_min, V_max, V_average;

  /**
   * @brief      Compute the min, max and average number of neighbor
   * Also compute the minimum distance in smoothing length
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_neighbors_stats(
      std::vector<body>& bodies,
      int64_t totalnbodies )
  {

    N_min = bodies[0].getNeighbors(); 
    uint64_t N_total = 0;
    for(auto& b: bodies)
    {
      uint64_t N = b.getNeighbors();
      N_min = std::min(N_min,N);
      N_max = std::max(N_max,N);
      N_total += N;
    }

    reduce_sum(N_total);
    reduce_min(N_min);
    reduce_max(N_max);
    reduce_min(min_dist);
    reduce_sum(average_dist_in_h);
    reduce_max(N_ghosts);
    average_dist_in_h /= totalnbodies;
    N_average = N_total/ totalnbodies;
  }

  /**
   * @brief      Compute the min, max and average smoothing length
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_smoothinglength_stats(
      std::vector<body>& bodies,
      int64_t totalnbodies)
  {
    double h_total = 0.;
    h_min = std::numeric_limits<double>::max();
    h_max = std::numeric_limits<double>::min();
    for(auto& b: bodies)
    {
      double h = b.radius();
      h_total += h;
      h_min = std::min(h,h_min);
      h_max = std::max(h,h_max);
    }
    reduce_sum(h_total);
    reduce_min(h_min);
    reduce_max(h_max);
    h_average = h_total / totalnbodies;
  }

  /**
   * @brief Compute the min, max and average norm of velocity
   */
  void
  compute_velocity_stats(
      std::vector<body>& bodies,
      int64_t totalnbodies)
  {
    V_min = std::numeric_limits<double>::max();
    V_max = std::numeric_limits<double>::min();
    V_average = 0.;
    double V_tot = 0.;
    for(auto& b: bodies)
    {
      double V = norm_point(b.getVelocity());
      V_max = std::max(V,V_max);
      V_min = std::min(V,V_min);
      V_tot += V;
    }
    V_average = V_tot/totalnbodies;
  }

  /**
   * @brief Periodic file output
   */
  void
  output(body_system<double,gdimension>& bs, const int rank)
  {
    static bool first_time = true;
    if (param::out_diagnostic_every <= 0
      || physics::iteration % param::out_diagnostic_every!=0)
      return;

    // compute diagnostic quantities
    bs.get_all(compute_neighbors_stats,bs.getNBodies());
    bs.get_all(compute_smoothinglength_stats,bs.getNBodies());
    bs.get_all(compute_velocity_stats,bs.getNBodies());

    // output only from rank #0
    if (rank != 0) return;
    const char * filename = "diagnostic.dat";

    if (first_time) {
      // Generate and output the header
      std::ostringstream oss_header;
      oss_header
        << "# Diagnostic: " <<std::endl
        << "# 1:iteration 2:time 3:h_min 4:h_max 5:h_avg "
        <<"6:N_min 7:N_max 8:N_avg 9:min_dist 10:avg_dist_h "
        <<"11:V_min 12:V_max 13:V_avg 14:N_ghosts" << std::endl;

      std::ofstream out(filename);
      out << oss_header.str();
      out.close();

      first_time = false;
    }

    std::ostringstream oss_data;
    oss_data << std::setw(5) << physics::iteration
      << std::setw(20) << std::scientific << std::setprecision(12)
      << physics::totaltime << std::setw(20)
      << h_min << std::setw(20) << h_max << std::setw(20)
      << h_average << std::setw(5) << N_min << std::setw(5)
      << N_max << std::setw(5) << N_average << std::setw(20)
      << min_dist << std::setw(20) << average_dist_in_h << std::setw(20)
      << V_min << std::setw(20) << V_max << std::setw(20) << V_average
      << std::setw(20) << N_ghosts;
    oss_data << std::endl;

    // Open file in append mode
    std::ofstream out(filename,std::ios_base::app);
    out << oss_data.str();
    out.close();

  } // scalar output

}; // namespace diagnostic

#endif // _PHYSICS_DIAGNOSTIC_H_
