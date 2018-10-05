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

// Define the relaxation 
// 1 = relaxation = computation of rot force 
// 0 = non relaxation = rotation applied
#define RELAXATION 1
#if RELAXATION == 1
#warning CAUTION RELAXATION MODE
#endif
#define ADIABATIC

#include "params.h"
#include "bodies_system.h"
#include "default_physics.h"
#include "BNS_physics.h"
#include "analysis.h"
#include "eos.h"
//#include "star_tracker.h"

namespace flecsi{
namespace execution{

static std::string initial_data_file;  // = initial_data_prefix  + ".h5part"
static std::string output_h5data_file; // = output_h5data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  // set kernel
  kernels::select(sph_kernel);

  // filenames (this will change for multiple files output)
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();
  oss << output_h5data_prefix << ".h5part";
  output_h5data_file = oss.str();

  // iteration and time
  physics::iteration = initial_iteration;
  physics::totaltime = initial_time;
  physics::dt = initial_dt; // TODO: use particle separation and Courant factor

  // set equation of state
  eos::select(eos_type);
}

void mpi_init_task(const char * parameter_file){
  using namespace param;
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  // set simulation parameters
  param::mpi_read_params(parameter_file);
  set_derived_params();

  // remove output file
  remove(output_h5data_file.c_str());

  // read input file
  body_system<double,gdimension> bs;
  bs.read_bodies(initial_data_file.c_str(),initial_iteration);
  bs.update_iteration();

  // auxiliary variables to record elapsed OMP time
  double wt = 0;
  double wt_start = 0;

  // set additional physics quantites
  physics::A = 0.6366197723675814;
  physics::angular_moment = bs.get_attribute<double>(
      initial_data_file.c_str(), "angularMomentum");

  bs.setMacangle(fmm_macangle);
  bs.setMaxmasscell(fmm_max_cell_mass);

  point_t momentum = {};

  if(physics::iteration == 0){
    // Set adiabatic value 
    bs.apply_all([](body_holder* source){
        source->getBody()->setAdiabatic(physics::A);
        physics::compute_internal_energy_from_adiabatic(source);
      });

  }else{
    rank|| clog(info) << "Recomputing A from u"<<std::endl; 
    rank|| clog(info) << "Considering velocityHalf = velocity" << std::endl;
    // Convert internal energy to A ratio 
    bs.apply_all(physics::compute_adiabatic_from_internal_energy);
    bs.apply_all([](body_holder* source){
        source->getBody()->setVelocityhalf(
          source->getBody()->getVelocity());
      });
  }

#ifdef OUTPUT
  if(out_scalar_every > 0 && physics::iteration % out_scalar_every == 0){
    MPI_Barrier(MPI_COMM_WORLD);
    bs.update_iteration();
    bs.apply_in_smoothinglength(physics::compute_density_pressure_soundspeed);
    
    // Compute conserved quantities
    bs.get_all(analysis::compute_lin_momentum);
    bs.get_all(analysis::compute_total_mass);
    bs.get_all(analysis::compute_total_energy);
    bs.get_all(analysis::compute_total_ang_mom);
    analysis::scalar_output("scalar_reductions.dat");
  }

  bs.write_bodies(output_h5data_prefix,physics::iteration);
#endif

  ++(physics::iteration); 
  do
  { 
    wt_start = omp_get_wtime();
    analysis::screen_output(rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // Integration step
    // TODO: leapfrog integration implemented incorrectly;
    //       update similarly to what is in the drivers/hydro.
    rank|| clog(trace)<<"Leapfrog integration"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.apply_all(physics::leapfrog_integration);
    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;

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
    rank|| clog(trace)<<"Rotation"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.apply_all(physics::apply_rotation);
    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;
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
    rank|| clog(trace)<<"compute_density_pressure_soundspeed"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.apply_in_smoothinglength(
      physics::compute_density_pressure_soundspeed); 
    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;

    bs.update_neighbors();

    // reset acceleration to 0 
    bs.apply_all([](body_holder* source){
      source->getBody()->setAcceleration(point_t{});
    });


    rank|| clog(trace)<<"Accel FMM"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.gravitation_fmm();
    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;

    rank|| clog(trace)<<"Accel hydro"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;

    rank|| clog(trace)<<"dadt"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.apply_in_smoothinglength(physics::compute_dadt);
    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;

#if RELAXATION == 1
    rank|| clog(trace)<<"Accel rot"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.get_all(physics::compute_QZZ);
    // Reduce QZZ
    MPI_Allreduce(MPI_IN_PLACE,&physics::QZZ,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    bs.get_all(physics::compute_rotation);
    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;
 #endif

    // Compute the new DT 
    physics::dt = 1.0;
    rank|| clog(trace)<<"DT computation"<<std::flush; 
    wt = omp_get_wtime(); 
    bs.apply_in_smoothinglength(physics::compute_dt);
    mpi_utils::reduce_min(physics::dt);

    rank|| clog(trace)<<".done "<< omp_get_wtime() - wt << "s" <<std::endl;
    rank|| clog(trace)<<"dt="<<physics::dt<<std::endl;
    assert(physics::dt != 1.0);

    // Compute internal energy for output 
    bs.apply_all(physics::compute_internal_energy_from_adiabatic);
    rank|| clog(trace) << "Total time=" << physics::totaltime << "s / "
                    << final_time << "s" << std::endl;

    if(out_scalar_every > 0 && physics::iteration % out_scalar_every == 0){
      // Compute the analysis values based on physics
      bs.get_all(analysis::compute_lin_momentum);
      bs.get_all(analysis::compute_total_mass);
      bs.get_all(analysis::compute_total_energy);
      bs.get_all(analysis::compute_total_ang_mom);
      // Only add the header in the first iteration
      analysis::scalar_output("scalar_reductions.dat");
    }

#ifdef OUTPUT
    if(out_h5data_every > 0 && physics::iteration % out_h5data_every == 0){
      bs.write_bodies(output_h5data_prefix,physics::iteration/out_h5data_every);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    ++(physics::iteration);
    physics::totaltime += physics::dt;
    
    MPI_Barrier(MPI_COMM_WORLD);
    rank|| clog(trace)<<"Iteration time = "<<omp_get_wtime()-wt_start<< "s" 
      <<std::endl;
  } while(physics::iteration < final_iteration);
}

flecsi_register_mpi_task(mpi_init_task,flecsi::execution);

void usage(int rank) {
  rank|| clog(warn)<<"Usage:"<<std::endl<<
    "./bns_3D [Starting iteration] [Max interation] [MAC Angle]"<<std::endl;
}

void
specialization_tlt_init(int argc, char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  rank|| clog(trace) << "In user specialization_driver" << std::endl;

  // check options list: exactly one option is allowed
  if (argc != 2) {
    rank|| clog(error) << "ERROR: parameter file not specified!" << std::endl;
    usage(rank);
    return;
  }

  flecsi_execute_mpi_task(mpi_init_task, flecsi::execution, argv[1]);

} // specialization driver

void 
driver(int argc,  char * argv[]){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  rank|| clog(warn) << "In user driver" << std::endl;
} // driver


} // namespace execution
} // namespace flecsi


