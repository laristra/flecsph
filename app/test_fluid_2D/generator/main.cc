/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 #include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "hdf5ParticleIO.h"


const double m_ = 1.0;
const double rho_ = 0.0;
const double u_ = 2.0; 
const double smoothing_length = 1e-2;
const char* fileprefix = "hdf5_fluid_2D";
const double timestep = 1e-3;
int32_t dimension = 2;

int main(int argc, char * argv[]){

  MPI_Init(&argc,&argv);
  
  int npart = 100;


  if(argc!=2){
    printf("./fluid_generator nparticles\n");
  }else{
    npart = atoi(argv[1]);
  }

  // Position
  double* x = new double[npart]();
  double* y = new double[npart]();
  // Velocity
  double* vx = new double[npart]();
  double* vy = new double[npart]();
  // Acceleration
  double* ax = new double[npart]();
  double* ay = new double[npart]();
  // Smoothing length 
  double* h = new double[npart]();
  // Density 
  double* rho = new double[npart]();
  // Pressure
  double* P = new double[npart]();
  // Mass
  double* m = new double[npart]();
  // Id
  int64_t* id = new int64_t[npart]();
  // Timestep 
  double* dt = new double[npart]();
  
  srand(time(NULL));

  std::cout<<"Generating: "<<npart<<" particles"<<std::endl;

  for(int i = 0 ; i < npart ; ++i){
    x[i] = (double)rand()/(double)RAND_MAX;
    y[i] = (double)rand()/(double)RAND_MAX;
    m[i] = m_;
    h[i] = smoothing_length;
    id[i] = i;
  }

  char filename[128];
  sprintf(filename,"%s.h5part",fileprefix);
  remove(filename);

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet; 
  testDataSet.createDataset(filename,MPI_COMM_WORLD);
    
  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",npart);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",dimension);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  //char * simName = "fluid_2D";
  //testDataSet.writeDatasetAttributeArray("name","string",simName);
  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",npart,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",npart,y);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",npart,vx);
  _d2.createVariable("vy",Flecsi_Sim_IO::point,"double",npart,vy);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("ax",Flecsi_Sim_IO::point,"double",npart,ax);
  _d2.createVariable("ay",Flecsi_Sim_IO::point,"double",npart,ay);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();


  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",npart,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",npart,rho);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",npart,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",npart,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",npart,id);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  testDataSet.closeFile(); 
  
  delete[] x;
  delete[] y;
  delete[] vx;
  delete[] vy;
  delete[] ax;
  delete[] ay;
  delete[] h;
  delete[] rho;
  delete[] P;
  delete[] m;
  delete[] id;
  delete[] dt;

  MPI_Finalize();

  return 0;
}
