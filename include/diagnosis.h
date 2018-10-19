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
 * @file diagnosis.h
 * @author Julien Loiseau
 * @date October 2018
 * @brief Basic diagnosis functions
 */

#ifndef _PHYSICS_DIAGNOSIS_H_
#define _PHYSICS_DIAGNOSIS_H_

#include <vector>
#include "params.h"

namespace diagnosis{

  uint64_t N_min, N_max, N_average; 
  double h_min, h_max, h_average; 
  double min_dist, average_dist_in_h;
  entity_id_t id_N_min, id_N_max, id_h_min, id_h_max;
  double V_min, V_max, V_average;

  /**
   * @brief      Compute the min, max and average number of neighbor
   * Also compute the minimum distance in smoothing length 
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_neighbors(
      std::vector<body_holder*>& bodies, tree_topology_t* tree, 
      int64_t totalnbodies )
  {
    min_dist = std::numeric_limits<double>::max(); 
    int64_t ncritical = 32; 
    uint64_t N_total = 0;
    average_dist_in_h = 0.;
    N_min = std::numeric_limits<std::uint64_t>::max();
    N_max = 0;
    // Find in radius the number of neighbors 
    tree->apply_sub_cells(
        tree->root(),
        0.,
        ncritical,
        param::sph_variable_h, 
        [](body_holder * srch, 
          std::vector<body_holder*>& nbh,
          double& min_d, double& average_dist)
        {
            // \TODO assert all bodies not nullptr and in interaction? 
            srch->getBody()->set_neighbors(nbh.size()-1); 
            
            double dist = std::numeric_limits<double>::max();
            double total_distance = 0.;
            for(auto n: nbh)
            {
              double d = distance(srch->coordinates(),
                n->coordinates());

              total_distance+=d; 

              if(d!=0.)
                dist = std::min(dist,d);
            }
            // Min distance 
            #pragma omp critical
            {
              min_d = std::min(min_d,dist);
              average_dist += total_distance/(nbh.size()-1);
            } 
        },min_dist,average_dist_in_h
    );
    for(auto b: bodies)
    {
      uint64_t N = b->getBody()->neighbors(); 
      N_min = std::min(N_min,N); 
      N_max = std::max(N_max,N);

      N_total += N; 
    }
    reduce_sum(N_total);
    reduce_min(N_min);
    reduce_max(N_max);
    reduce_min(min_dist);
    reduce_sum(average_dist_in_h);
    average_dist_in_h /= totalnbodies; 
    N_average = N_total/ totalnbodies;
  }

  /**
   * @brief      Compute the min, max and average smoothing length 
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_smoothinglength(
      std::vector<body_holder*>& bodies, int64_t totalnbodies)
  {
    double h_total = 0.;
    h_min = std::numeric_limits<double>::max();
    h_max = std::numeric_limits<double>::min();
    for(auto& b: bodies)
    {
      double h = b->getBody()->getSmoothinglength();
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
  compute_velocity(
      std::vector<body_holder*>& bodies, int64_t totalnbodies)
  {
    V_min = std::numeric_limits<double>::max();
    V_max = std::numeric_limits<double>::min();
    V_average = 0.;
    double V_tot = 0.;
    for(auto b: bodies)
    {
      double V = norm_point(b->getBody()->getVelocity());
      V_max = std::max(V,V_max);
      V_min = std::min(V,V_min);
      V_tot += V;
    }
    V_average = V_tot/totalnbodies;
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
    rank || clog(info)<< "#-- iteration:               time:" <<std::endl;
    rank || clog(info)
        << std::setw(14) << physics::iteration
        << std::setw(20) << std::scientific << std::setprecision(12)
        << physics::totaltime << std::endl;
    rank || clog(info)<<" h_min h_max h_avg N_min N_max N_avg"<< std::endl;
    rank || clog(info)<<h_min<<" "<<h_max<<" "<<h_average<<" "
      <<N_min<<" "<<N_max<<" "<<N_average<<
      std::endl;
  } 

  /**
   * @brief Periodic file output
   */
  void
  output(const char * filename = NULL)
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
        << "# Diagnosis: " <<std::endl
        << "# 1:iteration 2:time 3:h_min 4:h_max 5:h_avg "
        <<"6:N_min 7:N_max 8:N_avg 9:min_dist 10:avg_dist_h "
        <<"11:V_min 12:V_max 13:V_avg" << std::endl;

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
      << V_min << std::setw(20) << V_max << std::setw(20) << V_average;
    oss_data << std::endl;

    // Open file in append mode
    std::ofstream out(filename,std::ios_base::app);
    out << oss_data.str();
    out.close();

  } // scalar output

}; // diagnosis

#endif // _PHYSICS_DIAGNOSIS_H_
