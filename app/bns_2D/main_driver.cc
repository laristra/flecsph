/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*\

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

namespace flecsi{
namespace execution{

void
mpi_init_task(int startiteration, double macangle = 0){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int iter = startiteration; 
  int noutput = startiteration+1;
  int maxiter = 500;

  body_system<double,gdimension> bs;
  double maxtime = 1.0;
  double outputtime = 0.02;
  
  // Use the user reader defined in the physics file 
  bs.read_bodies("hdf5_bns_2D.h5part",startiteration);

  remove("output_bns_2D.h5part");
  const char outputname[64] = "output_bns_2D";

  bs.update_iteration();
    // For OMP time, just used on 0, initialized for others
  double start = 0;
  double start_iteration = 0;
  physics::dt = 0.001;
  bs.setMacangle(macangle);

  std::cout<<"MacAngle="<<macangle<<std::endl;

  point_t momentum = {};

  // Compute density, pressure, cs for next iteration
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

    // Compute and prepare the tree for this iteration 
    // - Compute the Max smoothing length 
    // - Compute the range of the system using the smoothinglength
    // - Cmopute the keys 
    // - Distributed qsort and sharing 
    // - Generate and feed the tree
    // - Exchange branches for smoothing length 
    // - Compute and exchange ghosts in real smoothing length 
    bs.update_iteration();

    // reset acceleration to 0 
    bs.apply_all([](body_holder* source){
      source->getBody()->setAcceleration(point_t{});
    });

    // Compute the hydro force 


    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Accel hydro"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Accel hydro"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_in_smoothinglength(physics::compute_internalenergy);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
    
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
      std::cout<<"Internal energy integration"<<std::flush; 
      start = omp_get_wtime(); 
    }
      bs.apply_all(physics::dudt_integration);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
    

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Leapfrog integration"<<std::flush; 
      start = omp_get_wtime(); 
    }
    if(iter == 1){
      bs.apply_all(physics::leapfrog_integration_first_step);
    }else{
      bs.apply_all(physics::leapfrog_integration);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;

    // Pressure density and soundspeed
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

    
    // Compute momentum 
    point_t total_momentum;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Momentum computation"<<std::flush; 
      start = omp_get_wtime(); 
    }
      bs.get_all(physics::compute_lin_momentum,&total_momentum);
    MPI_Barrier(MPI_COMM_WORLD);
    // Gather on 0 and sum 
    std::vector<point_t> sub_momentum;
    if(rank==0){
      sub_momentum.resize(size);
    }
    MPI_Gather(&total_momentum,sizeof(point_t),MPI_BYTE,
         &sub_momentum[0],sizeof(point_t),MPI_BYTE,0,MPI_COMM_WORLD);
    // Sum 
    if(rank==0){
      total_momentum = {0};
      for(auto v: sub_momentum){
        total_momentum += v;
      }
      // Display diff 
      std::cout << "="<< magnitude(point_to_vector(total_momentum-momentum)) << " ";

      momentum = total_momentum;
    }
    if(rank==0)
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
 
    // Compute the new DT 
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<"Accel hydro"<<std::flush; 
      start = omp_get_wtime(); 
    }
    bs.apply_in_smoothinglength(physics::compute_dt);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<".done "<< omp_get_wtime() - start << "s" <<std::endl;
      std::cout<<"dt="<<physics::dt<<std::endl;
    }

    physics::totaltime += physics::dt;

    if(rank == 0){
      std::cout<<"Total time="<<physics::totaltime<<"s / "<<maxtime<<"s"<<std::endl;
    }
#ifdef OUTPUT
    if(physics::totaltime > noutput*outputtime){
      bs.write_bodies(outputname,noutput);
      noutput++;
    }

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
  double macangle = 0;
  if(argc == 2){
    startiteration = atoi(argv[1]);
  }
  if(argc == 3){
    startiteration = atoi(argv[1]);
    macangle = atof(argv[2]);
  }

  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/
  /*std::string  filename(argv[1]);
  std::cout<<filename<<std::endl;*/
  flecsi_execute_mpi_task(mpi_init_task,startiteration,macangle); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


