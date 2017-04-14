#ifndef tree_driver_h
#define tree_driver_h

#include <iostream>

#include "physics.h"
#include "io.h"

namespace flecsi{
namespace execution{

int32_t serial_runtime_driver(int argc, char * argv[]){
  std::cout << "Serial runtime driver" << std::endl;
  
  if(argc != 3){
    std::cerr << "Usage: ./BNS_flecsi <dataFile> <nSteps>" << std::endl;
    exit(EXIT_FAILURE);
  }
  const char * filename = argv[1];
  const size_t nsteps = atoi(argv[2]);

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
   
    std::cout << "min: " << minposition << " max: " << maxposition << std::endl; 
    t.update_all(minposition,maxposition);
    
    // Compute density 
    for(auto bi : bodies){
      auto ents = t.find_in_radius(pool, bi->coordinates(), 2*bi->getSmoothinglength());
      std::vector<body*> bod = ents.to_vec();
      physics::computeDensity(bi,bod);
    } // for
   
    // Compute pressure 
    for(auto bi : bodies){
      physics::computePressure(bi);
      physics::computeSoundspeed(bi); 
    }
  
    // Compute Hydro 
    for(auto bi : bodies){
      auto ents = t.find_in_radius(pool, bi->coordinates(), 2*bi->getSmoothinglength());
      std::vector<body*> bod = ents.to_vec();
      physics::computeHydro(bi,bod);
    } // for

    // Compute grav
    for(auto bi: bodies){
      physics::computeGrav(bi,bodies);
    }

    // Compute acceleration 
    for(auto bi : bodies){
      physics::computeAcceleration(bi,dt);
    } 

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

    // Display first and last body
    std::cout << *bodies.front() << std::endl;
    std::cout << *bodies.back() << std::endl; 
  }
  return 0; 
} // driver

void
specialization_driver(
  int argc, 
  char * argv[]
)
{
  std::cout << "MPI runtime specialization_driver" << std::endl;
  //return 0;
} // specialization_driver

void 
driver(int argc, char * argv[]){
  std::cout << "MPI runtime driver" << std::endl;
}


} // namespace
} // namespace


#endif // tree_driver_h
