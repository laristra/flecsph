#include "io.h"

namespace io{

// Read data from txt file 
// This is not the most beautiful way, but enough for testing
void inputDataTxt(
    std::vector<body*>& bodies, 
    const char * filename,
    tree_topology_t& t
  ){
  // Load data from text file and setup bodies 
  FILE * inputfile = fopen(filename,"r");
  size_t nbodies = 0;
  double x, y, z, velX, velY, velZ, velHX, velHY, velHZ;
  double accX, accY, accZ, smoothinglength, pressure;
  double entropy, density, mass, tmp, angularMoment;
  int tmpD;
  while(fscanf(inputfile,
      "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
      " %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf \n",
      &x,&y,&z,                       // positions
      &velX,&velY,&velZ,              // velocity
      &velHX,&velHY,&velHZ,           // velHalf 
      &accX,&accY,&accZ,              // acceleration
      &density,&pressure,&entropy,    // density, pressure, entropy 
      &mass,&smoothinglength,&tmp,    // mass, smoothing length, dt 
      &tmpD,&tmp,&tmp,                // step, totalTime, tmp
      &tmp,&tmp,&angularMoment        // tmp, tmp, angularMoment
  )!=EOF){
    point_t position = {x,y,z}; 
    point_t velocity = {velX,velY,velZ};
    point_t velocityhalf = {velHX,velHY,velHZ};
    point_t acceleration = {accX, accY, accZ};
    auto bi = t.make_entity(position,velocity,velocityhalf,
        acceleration,density,pressure,entropy,mass,smoothinglength);   
    bodies.push_back(bi);
    ++nbodies;
    if(nbodies < 5)
      std::cout << *bi << std::endl;
  } // while 
} // inputData 

} // io
