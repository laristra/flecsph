#ifndef tree_driver_h
#define tree_driver_h

#include <iostream>
#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

#include "physics.h"
#include "io.h"

namespace flecsi{
namespace execution{

void
mpi_task()
{
  //int rank = 0;
  //int size = 0;
  //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //MPI_Comm_size(MPI_COMM_WORLD,&size);
  //std::cout << "Salut" <<  std::endl;
  //return rank;
} // mpi_task

flecsi_register_task(mpi_task,mpi,single);


//int32_t serial_runtime_driver(int argc, char * argv[]){
void specialization_driver(int argc, char * argv[]){
  std::cout << "Serial runtime driver" << std::endl;
  
  if(argc != 3){
    std::cerr << "Usage: ./BNS_flecsi <dataFile> <nSteps>" << std::endl;
    exit(EXIT_FAILURE);
  }
  const char * filename = argv[1];
  const size_t nsteps = atoi(argv[2]);

  auto dc = new flecsi::data::data_client_t;

  flecsi_register_data(*dc,hydro,pressure,double,global,1);
  flecsi_register_data(*dc,hydro,density,double,global,1);
  flecsi_register_data(*dc,hydro,position,point_t,global,1);
  flecsi_register_data(*dc,hydro,smoothinglength,double,global,1);

  auto accs = flecsi_get_accessors(*dc,hydro,double,global,0);

  //std::cout << "Launching task"; 
  //auto res = flecsi_execute_task(mpi_task,mpi,single);
  //std::cout << ".OK" << std::endl;
  //std::cout << "Wait" << std::endl;
  //res.wait();
  //std::cout << ".OK" << std::endl; 

  double timers[10] = {0.};
  tree_topology_t t; 
  std::vector<body*> bodies; 
  point_t maxposition = {-999,-999,-999}; 
  point_t minposition = {999,999,999};
  // Thread pool for the neighbor searh
  flecsi::thread_pool pool; 
  pool.start(16);

  // Load data from text file and setup bodies 
  io::inputDataTxt(bodies,filename,t); 

  // Compute max and min 
  for(auto bi : bodies){
    for(size_t i=0;i<bi->getPosition().dimension;++i){
     if(minposition[i] > bi->getPosition()[i])
       minposition[i] = bi->getPosition()[i];
     if(maxposition[i] < bi->getPosition()[i])
       maxposition[i] = bi->getPosition()[i];
   } // for
  } // for

  std::cout << "min: " << minposition << " max: " << maxposition << std::endl; 
  t.update_all(minposition,maxposition);


  // Add the bodies in the tree 
  for(auto bi : bodies){
    t.insert(bi);
  } // while

  double dt = 1.0e-7;

  for(size_t step = 0; step < nsteps; ++step){
    std::cout << "---- step = " << step+1 << "/" << nsteps << " dt: " << dt << std::endl;
  

    double step_01 = omp_get_wtime();
    {
    minposition = {999.,999.,999.};
    maxposition = {-999.,-999.,-999.};
    // move particles and compute min and max positions
    for(auto bi: bodies){
      physics::moveBody(bi,dt);
      // TODO add rotation
      for(size_t i=0;i<bi->getPosition().dimension;++i){
        if(minposition[i] > bi->getPosition()[i])
          minposition[i] = bi->getPosition()[i];
        if(maxposition[i] < bi->getPosition()[i])
          maxposition[i] = bi->getPosition()[i];
      } // for
      t.update(bi);
    } // for
    }
    double step_02 = omp_get_wtime();
    timers[0] += step_02 - step_01; 
    std::cout << "moveBody         : " << step_02-step_01 << "s" << std::endl; 
    {
    std::cout << "min: " << minposition << " max: " << maxposition << std::endl; 
    t.update_all(minposition,maxposition);
    }
    double step_03 = omp_get_wtime();
    timers[1] += step_03 - step_02;
    std::cout << "updateTree       : " << step_03-step_02 << "s" << std::endl; 
    {
    // Compute density 
    for(auto bi : bodies){
      auto ents = t.find_in_radius(pool, bi->coordinates(), 2*bi->getSmoothinglength());
      std::vector<body*> bod = ents.to_vec();
      physics::computeDensity(bi,bod);
    } // for
    }
    double step_04 = omp_get_wtime();
    timers[2] += step_04 - step_03;
    std::cout << "computeDensity   : " << step_04-step_03 << "s" << std::endl; 
    {
    // Compute pressure 
    for(auto bi : bodies){
      physics::computePressure(bi);
      physics::computeSoundspeed(bi); 
    }
    }
    double step_05 = omp_get_wtime();
    timers[3] += step_05 - step_04;
    std::cout << "computePressure  : " << step_05-step_04 << "s" << std::endl; 
    {
    // Compute Hydro 
    for(auto bi : bodies){
      auto ents = t.find_in_radius(pool, bi->coordinates(), 2*bi->getSmoothinglength());
      std::vector<body*> bod = ents.to_vec();
      physics::computeHydro(bi,bod);
    } // for
    }
    double step_06 = omp_get_wtime();
    timers[4] += step_06 - step_05;
    std::cout << "computeHydro     : " << step_06-step_05 << "s" << std::endl; 
    {
    // Compute grav
    for(auto bi: bodies){
      physics::computeGrav(bi,bodies);
    }
    }
    double step_07 = omp_get_wtime();
    timers[5] += step_07 - step_06;
    std::cout << "computeGrav      : " << step_07-step_06 << "s" << std::endl; 
    {
    // Compute acceleration 
    for(auto bi : bodies){
      physics::computeAcceleration(bi,dt);
    } 
    }
    double step_08 = omp_get_wtime();
    timers[6] += step_08 - step_07;
    std::cout << "computeAccel     : " << step_08-step_07 << "s" << std::endl; 
    {
    // Compute dt
    dt = 1.; 
    for(auto bi : bodies){
      point_t accel = bi->getAcceleration();
      double accelNorm = 0.;
      for(size_t i=0;i<accel.dimension;++i)
        accelNorm+=accel[i]*accel[i];
      accelNorm = sqrt(accelNorm);
      double ldt = physics::kCoeffDt*pow(bi->getSmoothinglength()/accelNorm,1./2.); 
      if(ldt < dt)
        dt = ldt;
    } // for
    }
    double step_09 = omp_get_wtime();
    timers[7] += step_09 - step_08;
    std::cout << "computeDt        : " << step_09-step_08 << "s" << std::endl; 
    // Display first and last body
    std::cout << *bodies.front() << std::endl;
    std::cout << *bodies.back() << std::endl; 
    std::cout << std::endl;
  }

  // Display mean times 
  std::cout << "moveBody         :" << timers[0]/(double)nsteps  << std::endl;
  std::cout << "updateTree       :" << timers[1]/(double)nsteps  << std::endl;
  std::cout << "computeDensity   :" << timers[2]/(double)nsteps  << std::endl;
  std::cout << "computePressure  :" << timers[3]/(double)nsteps  << std::endl;
  std::cout << "computeHydro     :" << timers[4]/(double)nsteps  << std::endl;
  std::cout << "computeGrav      :" << timers[5]/(double)nsteps  << std::endl;
  std::cout << "computeAccel     :" << timers[6]/(double)nsteps  << std::endl;
  std::cout << "computeDt        :" << timers[7]/(double)nsteps  << std::endl;

 // return 0; 
} // driver

//void
//specialization_driver(
//  int argc, 
//  char * argv[]
//)
//{
//  std::cout << "MPI runtime specialization_driver" << std::endl;
//  //return 0;
//} // specialization_driver

void 
driver(int argc, char * argv[]){
  std::cout << "MPI runtime driver" << std::endl;
}


} // namespace
} // namespace


#endif // tree_driver_h
