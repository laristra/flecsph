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
  
  /**
   * @brief      Compute the min, max and average number of neighbor
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_neighbors(
      std::vector<body_holder*>& bodies, tree_topology_t* tree, 
      int64_t totalnbodies )
  {
    int64_t ncritical = 32; 
    uint64_t N_total = 0;
    N_min = std::numeric_limits<std::uint64_t>::max();
    N_max = 0;
    // Find in radius the number of neighbors 
    tree->apply_sub_cells(
        tree->root(),
        0.,
        ncritical,
        param::sph_variable_h, 
        [](body_holder * srch, 
          std::vector<body_holder*>& nbh)
        {
            // \TODO assert all bodies not nullptr and in interaction? 
            srch->getBody()->set_neighbors(nbh.size()); 
        }
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
        <<"6:N_min 7:N_max 8:N_avg "<<std::endl;

      std::ofstream out(filename);
      out << oss_header.str();
      out.close();

      first_time = false;
    }

    std::ostringstream oss_data;
    oss_data << std::setw(14) << physics::iteration
      << std::setw(20) << std::scientific << std::setprecision(12)
      << physics::totaltime << std::setw(20) 
      << h_min << std::setw(20) << h_max << std::setw(20) 
      << h_average << std::setw(10) << N_min << std::setw(10)
      << N_max << std::setw(10) << N_average ; 
    oss_data << std::endl;

    // Open file in append mode
    std::ofstream out(filename,std::ios_base::app);
    out << oss_data.str();
    out.close();

  } // scalar output

}; // diagnosis

#endif // _PHYSICS_DIAGNOSIS_H_
