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
 * @file analysis.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic analysis
 */

#ifndef _PHYSICS_ANALYSIS_H_
#define _PHYSICS_ANALYSIS_H_

#include <vector>
#include "params.h"

// OpenMP point reduction
#pragma omp declare reduction(add_point : point_t : omp_out += omp_in) \
initializer(omp_priv=point_t{})

//#include "physics.h"

namespace analysis{

  enum e_conservation: size_t
  { MASS = 0 , ENERGY = 1, MOMENTUM = 2, ANG_MOMENTUM = 3 };

  point_t linear_momentum;
  point_t total_ang_mom;
  double total_mass;
  double total_energy;
  double total_kinetic_energy;
  double total_internal_energy;
  double velocity_part;

  /**
   * @brief      Compute the linear momentum
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_momentum(
      std::vector<body>& bodies)
  {
    linear_momentum = {0};
    #pragma omp parallel for reduction(add_point:linear_momentum)
    for(size_t i = 0 ; i < bodies.size(); ++i){
      if(bodies[i].type() != NORMAL)  continue;
      linear_momentum += bodies[i].mass()*bodies[i].getVelocity();
    }
    mpi_utils::reduce_sum(linear_momentum);
  }

  /**
   * @brief      Compute the total mass
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_mass(
      std::vector<body>& bodies)
  {
    total_mass = 0.;
    #pragma omp parallel for reduction(+:total_mass)
    for(size_t i = 0 ; i < bodies.size(); ++i) {
      if(bodies[i].type() != NORMAL)  continue;
      total_mass += bodies[i].mass();
    }
    mpi_utils::reduce_sum(total_mass);
  }

  /**
   * @brief      Compute the total energy (internal + kinetic)
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_energy(
      std::vector<body>& bodies)
  {
    using namespace param;

    total_energy = 0.;
    if (thermokinetic_formulation) {
      #pragma omp parallel for reduction(+:total_energy)
      for(size_t i = 0 ; i < bodies.size(); ++i){
        if(bodies[i].type() != NORMAL)  continue;
        total_energy += bodies[i].mass()*bodies[i].getTotalenergy();
      }
    }
    else {
      #pragma omp parallel for reduction(+:total_energy)
      for(size_t i = 0 ; i < bodies.size(); ++i){
        if(bodies[i].type() != NORMAL)  continue;
        double m = bodies[i].mass(),
            eint = bodies[i].getInternalenergy();
        total_energy += m*eint;
        point_t v = bodies[i].getVelocity();
        double v2 = v[0]*v[0];
        for(unsigned short int k=1; k<gdimension; ++k)
          v2 += v[k]*v[k];
        total_energy += .5*m*v2;
      }
    }
    mpi_utils::reduce_sum(total_energy);
  }

  /**
   * @brief      Sum up total kinetic energy
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_kinetic_energy(
      std::vector<body>& bodies)
  {
    using namespace param;

    total_kinetic_energy = 0.;
    #pragma omp parallel for reduction(+:total_kinetic_energy)
    for(size_t i = 0 ; i < bodies.size(); ++i){
      if(bodies[i].type() != NORMAL)  continue;
      double m = bodies[i].mass();
      point_t v = bodies[i].getVelocity();
      double v2 = v[0]*v[0];
      for(unsigned short int k=1; k<gdimension; ++k)
        v2 += v[k]*v[k];
      total_kinetic_energy += .5*m*v2;
    }
    mpi_utils::reduce_sum(total_kinetic_energy);
  }

  /**
   * @brief      Sum up internal energy
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_internal_energy(
      std::vector<body>& bodies)
  {
    using namespace param;

    total_internal_energy = 0.;
    #pragma omp parallel for reduction(+:total_internal_energy)
    for(size_t i = 0 ; i < bodies.size(); ++i) {
      if(bodies[i].type() != NORMAL)  continue;
      total_internal_energy += bodies[i].mass() * bodies[i].getInternalenergy();
    }
    mpi_utils::reduce_sum(total_internal_energy);
  }

  /**
   * @brief      Compute total angular momentum
   *
   * @param      bodies  Vector of all the local bodies
   */
  void
  compute_total_ang_mom(
      std::vector<body>& bodies)
  {
    total_ang_mom = {0};
    if constexpr (gdimension == 2) {
      #pragma omp parallel for reduction(add_point:total_ang_mom)
      for(size_t i = 0 ; i < bodies.size(); ++i){
        if(bodies[i].type() != NORMAL)  continue;
        const double m = bodies[i].mass();
        const point_t v = bodies[i].getVelocity();
        const point_t r = bodies[i].coordinates();
        total_ang_mom[0] += m*(r[0]*v[1] - r[1]*v[0]);
      }
      mpi_utils::reduce_sum(total_ang_mom);

    }
    else if constexpr (gdimension == 3) {
      #pragma omp parallel for reduction(add_point:total_ang_mom)
      for(size_t i = 0 ; i < bodies.size(); ++i){
        if(bodies[i].type() != NORMAL)  continue;
        const double m = bodies[i].mass();
        const point_t v = bodies[i].getVelocity();
        const point_t r = bodies[i].coordinates();
        total_ang_mom[0] += m*(r[1]*v[2] - r[2]*v[1]);
        total_ang_mom[1] += m*(r[2]*v[0] - r[0]*v[2]);
        total_ang_mom[2] += m*(r[0]*v[1] - r[1]*v[0]);
      }
      mpi_utils::reduce_sum(total_ang_mom);
    }
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
      (++count-1)%screen_length ||
      clog_one(info)<< "#-- iteration:               time:" <<std::endl;
      clog_one(info)
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
  scalar_output(body_system<double,gdimension>& bs, const int rank)
  {
    static bool first_time = true;
    if(param::out_scalar_every <= 0 ||
       physics::iteration % param::out_scalar_every != 0)
       return;

    // compute reductions
    bs.get_all(compute_total_momentum);
    bs.get_all(compute_total_mass);
    bs.get_all(compute_total_energy);
    bs.get_all(compute_total_kinetic_energy);
    bs.get_all(compute_total_internal_energy);
    bs.get_all(compute_total_ang_mom);

    // output only from rank #0
    if (rank != 0) return;
    const char * filename = "scalar_reductions.dat";

    if (first_time) {
      // Generate and output the header
      std::ostringstream oss_header;
      switch(gdimension) {
      case 1:
        oss_header
          << "# Scalar reductions: " <<std::endl
          << "# 1:iteration 2:time 3:timestep 4:total_mass "
          << " 5:total_energy 6:kinetic_energy 7:internal_energy"
          << " 8:mom_x"<< std::endl;
        break;

      case 2:
        oss_header
          << "# Scalar reductions: " <<std::endl
          << "# 1:iteration 2:time 3:timestep 4:total_mass 5:total_energy"
          <<  " 6:kinetic_energy 7:internal_energy" <<std::endl
          << "# 8:mom_x 9:mom_y 10:and_mom_z"<< std::endl;
        break;

      case 3:
      default:
        oss_header
          << "# Scalar reductions: " <<std::endl
          << "# 1:iteration 2:time 3:timestep 4:total_mass 5:total_energy"
          <<  " 6:kinetic_energy 7:internal_energy "  <<std::endl
          << "# 8:mom_x 9:mom_y 10:mom_z "
          <<  "11:ang_mom_x 12:ang_mom_y 13:ang_mom_z"<< std::endl;
      }

      std::ofstream out(filename);
      out << oss_header.str();
      out.close();
      first_time = false;
    }

    std::ostringstream oss_data;
    oss_data << std::setw(14) << physics::iteration
      << std::setw(20) << std::scientific << std::setprecision(12)
      << physics::totaltime << std::setw(20) << physics::dt <<" "
      << total_mass <<" "<< total_energy <<" "
      << total_kinetic_energy<<" "<<total_internal_energy<<" ";
    for(unsigned short int k = 0 ; k < gdimension ; ++k)
      oss_data <<" "<< linear_momentum[k];

    if(gdimension==2)
      oss_data <<" "<< total_ang_mom[0];

    if(gdimension==3)
      for(unsigned short k = 0 ; k < gdimension ; ++k)
        oss_data <<" "<< total_ang_mom[k];

    oss_data << std::endl;

    // Open file in append mode
    std::ofstream out(filename,std::ios_base::app);
    out << oss_data.str();
    out.close();

  } // scalar output


  bool
  check_conservation(
    const std::vector<e_conservation>& check
  )
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    clog_one(info) << "Checking conservation of: ";
    for(auto c: check){
      switch(c){
        case MASS:
          clog_one(info) << " MASS ";
          break;
        case ENERGY:
          clog_one(info) << " ENERGY ";
          break;
        case MOMENTUM:
          clog_one(info) << " MOMENTUM ";
          break;
        case ANG_MOMENTUM:
          clog_one(info) << " ANG_MOMENTUM ";
          break;
        default:
          break;
      } // switch
    }
    // Open file and check
    double tol = 1.0e-16;
    size_t iteration;
    double time, timestep, mass, energy, momentum[gdimension], ang_mom[3];
    double base_energy, base_mass, base_ang_mom[3], base_momentum[gdimension];
    std::ifstream inFile;
    inFile.open("scalar_reductions.dat");
    if (!inFile) {
      std::cerr << "Unable to open file scalar_reductions.dat";
      exit(1);   // call system to stop
    }

    // Read the first line
    inFile >> iteration >> time >> timestep >> base_mass >> base_energy;
    for(size_t i = 0 ; i < gdimension ; ++i){
      inFile >> base_momentum[i];
    }
    if(gdimension == 2){
      inFile >> base_ang_mom[2];
    }
    if(gdimension == 3){
      inFile >> base_ang_mom[0] >> base_ang_mom[1] >> base_ang_mom[2];
    }

    while (inFile.good()) {
      inFile >> iteration >> time >> timestep >> mass >> energy;
      for(size_t i = 0 ; i < gdimension ; ++i){
        inFile >> momentum[i];
      }
      if(gdimension == 2){
        inFile >> ang_mom[2];
      }
      if(gdimension == 3){
        inFile >> ang_mom[0] >> ang_mom[1] >> ang_mom[2];
      }
      // Check the quantities
      // Check only the required ones
      for(auto c: check){
        switch(c){
          case MASS:
            if(!(abs(base_mass - mass) < tol)){
              std::cerr<<"Mass is not conserved"<<std::endl;
              return false;
            }
            break;
          case ENERGY:
            if(!(abs(base_energy - energy) < tol)){
              std::cerr<<"Energy is not conserved"<<std::endl;
              return false;
            }
            break;
          case MOMENTUM:
            if(!(abs(base_momentum - momentum) < tol)){
              std::cerr<<"Momentum is not conserved"<<std::endl;
              return false;
            }
            break;
          case ANG_MOMENTUM:
            if(!(abs(base_ang_mom - ang_mom) < tol)){
              std::cerr<<"Angular Momentum is not conserved"<<std::endl;
              return false;
            }
            break;
          default:
            break;
        } // switch
      } // for
    } // while
    // Close the file
    inFile.close();
    return true;
  } // conservation check

}; // physics

#endif // _PHYSICS_ANALYSIS_H_
