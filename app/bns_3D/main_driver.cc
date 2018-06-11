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

#include <bodies_system.h>
#include "physics.h"

#define FROT 0
#define ROT 1

namespace flecsi{
namespace execution{

void
mpi_init_task(int startiteration = 0, int maxiter = 1000, double macangle = 0){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int iter = startiteration; 
  int noutput = 0;//startiteration+1;

  body_system<double,gdimension> bs;
  double maxtime = 1.0;
  double outputtime = 0.05;

  
  // Use the user reader defined in the physics file 
  bs.read_bodies("hdf5_bns_3D.h5part",startiteration);

  remove("output_bns_3D.h5part");
  const char outputname[64] = "output_bns_3D";

  bs.update_iteration();
    // For OMP time, just used on 0, initialized for others
  double start = 0;
  double start_iteration = 0;
  physics::dt = 1.0e-8;
  bs.setMacangle(macangle);
  bs.setMaxmasscell(10e-5);

  std::cout<<"MacAngle="<<macangle<<std::endl;

  point_t momentum = {};



  
  if(iter == 0){
    // Set adiabatic value 
    bs.apply_all([](body_holder* source){
        source->getBody()->setAdiabatic(physics::A);
        physics::compute_internal_energy(source);
      });

  }else{
    std::cout << "Recomputing A from u"<<std::endl; 
    std::cout << "Considering velocityHalf = velocity" << std::endl;
    // Convert internal energy to A ratio 
    bs.apply_all(physics::compute_adiabatic_from_u);
    bs.apply_all([](body_holder* source){
        source->getBody()->setVelocityhalf(
          source->getBody()->getVelocity());
      });
  }

#ifdef OUTPUT
  bs.write_bodies(outputname,startiteration);
#endif

  ++iter; 
  do
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    start_iteration = omp_get_wtime();
    if(rank==0)
      std::cout<<std::endl<<"#### Iteration "<<iter<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // Integration step
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Leapfrog integration"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_all(physics::leapfrog_integration);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    // Get rid of distant particles for faster run 
    bs.apply_all([](body_holder* source){
        body * src = source->getBody();
        point_t pos = src->getPosition();
        for(size_t i = 0 ; i < gdimension ; ++i){
          if(pos[i] > 10){
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

#if ROT == 1
    // Rotation of the stars
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Rotation"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_all(physics::apply_rotation);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
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
    if(rank==0){
      std::cout<<"compute_density_pressure_soundspeed"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_in_smoothinglength(
      physics::compute_density_pressure_soundspeed); 
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    bs.update_neighbors();

    // reset acceleration to 0 
    bs.apply_all([](body_holder* source){
      source->getBody()->setAcceleration(point_t{});
    });


    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Accel FMM"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.gravitation_fmm();
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Accel hydro"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

#if FROT == 1
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Accel rot"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.get_all(physics::compute_QZZ);
    // Reduce QZZ
    MPI_Allreduce(MPI_IN_PLACE,&physics::QZZ,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    bs.get_all(physics::compute_rotation);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
 #endif

    // Compute the new DT 
    physics::dt = 1.0;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"DT computation"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_in_smoothinglength(physics::compute_dt);
    // Get the minimum using MPI
    MPI_Allreduce(MPI_IN_PLACE,&physics::dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
      std::cout<<"dt="<<physics::dt<<std::endl;
    }
    assert(physics::dt != 1.0);
    physics::totaltime += physics::dt;

    // Compute internal energy for output 
    bs.apply_all(physics::compute_internal_energy);

    if(rank == 0){
      std::cout<<"Total time="<<physics::totaltime<<"s / "<<maxtime<<"s"<<std::endl;
    }
#ifdef OUTPUT
    if(physics::totaltime >= noutput*outputtime){
      bs.write_bodies(outputname,noutput);
      noutput++;
    }


    //if(iter % 50 == 0){
    //  bs.write_bodies(outputname,noutput);
    //  noutput++;
    //}

    //if(iter % iteroutput == 0){ 
    //  bs.write_bodies("output_fluid",iter/iteroutput);
    //}
#endif
    ++iter;
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
      std::cout<<"Iteration time = "<<omp_get_wtime()-start_iteration<< "s" <<std::endl;
    }
  }while(iter<maxiter);
  //}while(physics::totaltime<physics::maxtime);
}

flecsi_register_mpi_task(mpi_init_task);

void 
specialization_tlt_init(int argc, char * argv[]){
  
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

  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/
  /*std::string  filename(argv[1]);
  std::cout<<filename<<std::endl;*/
  flecsi_execute_mpi_task(mpi_init_task,startiteration,maxiter,macangle); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


