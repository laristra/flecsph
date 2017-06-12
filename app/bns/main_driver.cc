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

#include "bodies_system.h"
#include "tree_colorer.h"
#include "physics.h"
#include "io.h"


namespace flecsi{
 namespace execution{

void
mpi_init_task(int startiteration){
  const char * filename = "../data/data_bns_4169.txt";
  //const char * filename = "../data/data_binary_rdy_16288.txt";
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int totaliters = 100000;
  int iteroutput = 10;
  double totaltime = 0.0;
  double maxtime = 100.0;
  int iter = startiteration;
  
  physics::K = 0.63662;  

  body_system<double,gdimension> bs;
  bs.setMaxmasscell(1.0e-5);
  bs.setMacangle(0.7);
  
  bs.read_bodies_txt(filename); 
  //io::inputDataTxtRange(rbodies,nbodies,totalnbodies,rank,size,filename); 

#ifdef OUTPUT
  bs.write_bodies("output_bns.h5part",iter); 
#endif

  ++iter; 
  do
  {  
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<std::endl<<"#### Iteration "<<iter<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    
    bs.update_iteration();
     
    // Compute Density
    if(rank == 0) 
      std::cout<<"Density"<<std::flush;
    bs.apply_in_smoothinglength(physics::compute_density);
    if(rank == 0)
      std::cout<<".done"<<std::endl<<std::flush;
  
    // Compute Pressure and SoundSpeed
    if(rank == 0)
      std::cout<<"Pressure"<<std::flush;
    bs.apply_all(physics::compute_pressure);
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    if(rank == 0)
      std::cout<<"Soundspeed"<<std::flush;
    bs.apply_all(physics::compute_soundspeed);
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    bs.update_neighbors();
    
    // Compute Hydro 
    if(rank==0)
      std::cout<<"Hydro acceleration"<<std::flush;
    bs.apply_in_smoothinglength(physics::compute_hydro_acceleration);
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;
     
    // Compute and exchange the aggregate  
    // Then compute gravitational acceleration
    if(rank==0)
      std::cout<<"Grav acceleration"<<std::flush;
    bs.gravitation_fmm();
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;


    // Compute internal energy 
    if(rank==0)
      std::cout<<"Hydro acceleration"<<std::flush;
    bs.apply_in_smoothinglength(physics::compute_internalenergy);
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    
    // Leapfrog part one 
    if(rank==0)
      std::cout<<"MoveParticles"<<std::flush;    
    bs.apply_all(physics::leapfrog_integration_1);
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    // Compute the new DT
    physics::dt = 1;
    if(rank==0)
      std::cout<<"Dt"<<std::flush;
    bs.apply_in_smoothinglength(physics::compute_dt);
    
    // Do a reduction for dt
    MPI_Allreduce(MPI_IN_PLACE,&physics::dt,1,MPI_DOUBLE,
        MPI_MIN,MPI_COMM_WORLD);
 
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;
    assert(physics::dt<1&&physics::dt>0);

    if(rank==0)
      std::cout<<"dt="<<physics::dt<<std::endl<<std::flush;

    // Move the particles 
    if(rank==0)
      std::cout<<"MoveParticles"<<std::flush;    
    bs.apply_all(physics::leapfrog_integration_2);
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    // Add rotation 


#ifdef OUTPUT
    if(iter % iteroutput == 0)
    {
      bs.write_bodies("output_bns.h5part",iter/iteroutput);
    }
#endif
    ++iter;
    totaltime += physics::dt; 
    if(rank==0)
      std::cout<<"Time: "<<totaltime<<std::endl;

  }while(totaltime<maxtime);
} 

flecsi_register_task(mpi_init_task,mpi,index);

void 
specialization_driver(int argc, char * argv[]){
 
  int startiteration = 0;

  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/
  /*std::string filename(argv[1]);
  std::cout<<filename<<std::endl;*/
  flecsi_execute_task(mpi_init_task,mpi,index,startiteration); 
} // specialization driver

void 
driver(int argc, char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


