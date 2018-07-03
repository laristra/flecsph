/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*\

 /*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@ @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@///// /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@      /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@ /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////  /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@      /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@      @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file main_driver.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Specialization and Main driver used in FleCSI. 
 * The Specialization Driver is normally used to register data and the main 
 * code is in the Driver.  
 */

#include <iostream>
#include <numeric> // For accumulate
#include <iostream>

#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

#include "params.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "BNS_physics.h"
#include "analysis.h"

// Define the relaxation 
// 1 = relaxation = computation of rot force 
// 0 = non relaxation = rotation applied
#define RELAXATION 1
#if RELAXATION == 1
#warning CAUTION RELAXATION MODE
#endif

namespace flecsi{
namespace execution{

void
mpi_init_task(int startiteration = 0, int maxiter = 1000, double macangle = 0){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  clog_set_output_rank(0);
  
  physics::iteration = startiteration; 
  int noutput = 0;//startiteration+1;

  body_system<double,gdimension> bs;
  double maxtime = 1.0;
  double outputtime = 0.05;

  double outputtime_analysis = 0.02;
  int noutput_analysis = 0;

  
  // Use the user reader defined in the physics file 
  bs.read_bodies("hdf5_bns_3D.h5part",startiteration);

  remove("output_bns_3D.h5part");
  const char outputname[64] = "output_bns_3D";

  bs.update_iteration();
    // For OMP time, just used on 0, initialized for others
  double start = 0;
  double start_iteration = 0;

// Setting physics data 
  physics::dt = 1.0e-8;
  //physics::gamma = 2.0;  // converted to a parameter (poly_gamma)
  physics::A = 0.6366197723675814;
  // Load angular momentum from data file 
  physics::angular_moment = bs.get_attribute<double>("hdf5_bns_3D.h5part", 
    "angularMomentum");
//

  bs.setMacangle(macangle);
  bs.setMaxmasscell(10e-5); // TODO: hardcoded value, turn to a parameter

  clog_one(trace)<<"MacAngle="<<macangle<<std::endl;

  point_t momentum = {};

  if(physics::iteration == 0){
    // Set adiabatic value 
    bs.apply_all([](body_holder* source){
        source->getBody()->setAdiabatic(physics::A);
        physics::compute_internal_energy_from_adiabatic(source);
      });

  }else{
    clog_one(info) << "Recomputing A from u"<<std::endl; 
    clog_one(info) << "Considering velocityHalf = velocity" << std::endl;
    // Convert internal energy to A ratio 
    bs.apply_all(physics::compute_adiabatic_from_internal_energy);
    bs.apply_all([](body_holder* source){
        source->getBody()->setVelocityhalf(
          source->getBody()->getVelocity());
      });
  }

#ifdef OUTPUT
  bs.write_bodies(outputname,startiteration);
#endif

  ++(physics::iteration); 
  do
  { 
    start_iteration = omp_get_wtime();
    analysis::screen_output();
    MPI_Barrier(MPI_COMM_WORLD);

    // Integration step
    clog_one(trace)<<"Leapfrog integration"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_all(physics::leapfrog_integration);
    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    // Get rid of distant particles for faster run 
    bs.apply_all([](body_holder* source){
        body * src = source->getBody();
        point_t pos = src->getPosition();
        for(size_t i = 0 ; i < gdimension ; ++i){
          if(pos[i] > 10){ // TODO: 10 is hardcoded, turn to a parameter
            pos[i] = 10;
            src->setVelocity(point_t{});
          }
          if(pos[i] < -10){
            pos[i] = -10;
            src->setVelocity(point_t{});
          }
          src->setPosition(pos);
        }
    });

#if RELAXATION == 0
    // Rotation of the stars
    clog_one(trace)<<"Rotation"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_all(physics::apply_rotation);
    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
#endif

    // Compute and prepare the tree for this iteration 
    // - Compute the Max smoothing length 
    // - Compute the range of the system using the smoothinglength
    // - Cmopute the keys 
    // - Distributed qsort and sharing 
    // - Generate and feed the tree
    // - Exchange branches for smoothing length 
    // - Compute and exchange ghosts in real smoothing length 
    bs.update_iteration();

    MPI_Barrier(MPI_COMM_WORLD);
    clog_one(trace)<<"compute_density_pressure_adiabatic_soundspeed"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_in_smoothinglength(
      physics::compute_density_pressure_adiabatic_soundspeed); 
    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    bs.update_neighbors();

    // reset acceleration to 0 
    bs.apply_all([](body_holder* source){
      source->getBody()->setAcceleration(point_t{});
    });


    clog_one(trace)<<"Accel FMM"<<std::flush; 
    start = omp_get_wtime(); 
    bs.gravitation_fmm();
    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    clog_one(trace)<<"Accel hydro"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    clog_one(trace)<<"dadt"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_in_smoothinglength(physics::compute_dadt);
    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

#if RELAXATION == 1
    clog_one(trace)<<"Accel rot"<<std::flush; 
    start = omp_get_wtime(); 
    bs.get_all(physics::compute_QZZ);
    // Reduce QZZ
    MPI_Allreduce(MPI_IN_PLACE,&physics::QZZ,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    bs.get_all(physics::compute_rotation);
    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
 #endif

    // Compute the new DT 
    physics::dt = 1.0;
    clog_one(trace)<<"DT computation"<<std::flush; 
    start = omp_get_wtime(); 
    bs.apply_in_smoothinglength(physics::compute_dt);
    mpi_utils::reduce_min(physics::dt);

    clog_one(trace)<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
    clog_one(trace)<<"dt="<<physics::dt<<std::endl;
    assert(physics::dt != 1.0);
    physics::totaltime += physics::dt;

    // Compute internal energy for output 
    bs.apply_all(physics::compute_internal_energy_from_adiabatic);
    clog_one(trace)<<"Total time="<<physics::totaltime<<"s / "<<maxtime<<"s"
      <<std::endl;

#ifdef OUTPUT_ANALYSIS
    if(physics::totaltime >= noutput_analysis*outputtime_analysis){
      // Compute the analysis values based on physics 
      bs.get_all(analysis::compute_lin_momentum);
      bs.get_all(analysis::compute_total_mass);
      // Output 
      noutput_analysis++;
      // Only add the header in the first iteration
      analysis::scalar_output("scalar_bns_3D.dat");
    }
#endif



    // OUTPUT step
#ifdef OUTPUT
    if(physics::totaltime >= noutput*outputtime){
      bs.write_bodies(outputname,noutput);
      noutput++;
    }
#endif
    ++(physics::iteration);
    
    MPI_Barrier(MPI_COMM_WORLD);
    clog_one(trace)<<"Iteration time = "<<omp_get_wtime()-start_iteration<< "s" 
      <<std::endl;
  }while(physics::iteration<maxiter);
  //}while(physics::totaltime<physics::maxtime);
}

flecsi_register_mpi_task(mpi_init_task);

void 
usage()
{
  clog_one(warn)<<"Usage:"<<std::endl<<
    "./bns_3D [Starting iteration] [Max interation] [MAC Angle]"<<std::endl;
}

void 
specialization_tlt_init(int argc, char * argv[]){
  
  usage();

  // Default start at iteration 0
  int startiteration = 0;
  int maxiter = 1000;
  double macangle = 0;
  if(argc == 2){
    startiteration = atoi(argv[1]);
  }
  if(argc == 3){
    startiteration = atoi(argv[1]);
    maxiter = atoi(argv[2]);
  }
  if(argc == 4){
    startiteration = atoi(argv[1]);
    maxiter = atoi(argv[2]);
    macangle = atof(argv[3]);
  }

  clog_one(warn) << "In user specialization_driver" << std::endl;
  flecsi_execute_mpi_task(mpi_init_task,startiteration,maxiter,macangle); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  clog_one(warn) << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


