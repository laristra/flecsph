/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 #include <iostream>
#include <algorithm>
#include <cassert>

#include "hdf5ParticleIO.h"
#include "physics/kernel.h"


const double ldistance = 0.001;  // Distance between the particles 
const double localgamma = 5./3.;
const double rho_1 = 1;
//const double rho_2 = 0.125;
const double pressure_1 = 10e-7;
//const double pressure_2 = 0.1;
const double u_1 = .5;
//const double u_2 = 2;
const double m_1 = 1.0e-7;
//const double m_2 = 1.0e-5;
const double smoothing_length = 2*ldistance;
const char* fileprefix = "hdf5_sedov";

bool 
in_radius(
    double x,
    double y,
    double x0,
    double y0, 
    double r)
{
  return (x-x0)*(x-x0)+(y-y0)*(y-y0)<r*r;
}

double 
density(
    double x,
    double y,
    double x0,
    double y0)
{
  double radius = sqrt(pow(x-x0,2.)+pow(y-y0,2.)); 
  if(radius==0.){
    return 100.;
  }
  return (100.-radius)/(radius); 
}

int main(int argc, char * argv[]){


  int64_t sparticles = 100;
  if(argc!=2){
    printf("./sedov_generator [square nParticles]\n");
    fprintf(stderr,"Generating default number of particles=%ld*%ld=%ld",
        sparticles,sparticles,sparticles*sparticles);
  }else{
    sparticles = atoll(argv[1]);
  }

  int rank, size; 
  int provided; 
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE); 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int64_t nparticles = sparticles*sparticles;

  double radius = (ldistance*(sparticles-1))/2.; 
  double x_c = (sparticles-1)*ldistance/2.;
  double y_c = x_c;//(sparticles-1)*ldistance/2.;

  std::cout<<"Sphere: r="<<radius<<" pos=["<<x_c<<";"<<y_c<<"]"<<std::endl;

  //if(rank==0){
    printf("Generating %ld particles by %ldx%ld in sphere r=%.4f\n",
        nparticles,sparticles,sparticles,radius);
  //}

  // Start on  0 0
  double x_topproc = x_c-radius;
  double y_topproc = y_c-radius;

  double maxxposition = /*(sparticles-1)*ldistance;*/x_c+radius;
  double maxyposition = /*(sparticles-1)*ldistance;*/y_c+radius;
  
  // Position
  double* x = new double[nparticles]();
  double* y = new double[nparticles]();
  double* z = new double[nparticles]();
  // Velocity
  double* vx = new double[nparticles]();
  double* vy = new double[nparticles]();
  double* vz = new double[nparticles]();
  // Acceleration
  double* ax = new double[nparticles]();
  double* ay = new double[nparticles]();
  double* az = new double[nparticles]();
  // Smoothing length 
  double* h = new double[nparticles]();
  // Density 
  double* rho = new double[nparticles]();
  // Internal Energy 
  double* u = new double[nparticles]();
  // Pressure
  double* P = new double[nparticles]();
  // Mass
  double* m = new double[nparticles]();
  // Id
  int64_t* id = new int64_t[nparticles]();
  // Timestep 
  double* dt = new double[nparticles]();
  
  // Id of my first particle 
  int64_t posid = 0;

  // Header data 
  // the number of particles = nparticles 
  // The value for constant timestep 
  double timestep = 0.001;
  int dimension = 2;
  
  std::cout<<"top_X="<<x_topproc<<" top_Y="<<y_topproc<<
    " maxX="<<maxxposition<<" maxY="<<maxyposition<<std::endl;

  double xposition = /*0;*/x_topproc; 
  int64_t tparticles = 0;
  double yposition = /*0;*/y_topproc;
  //int xpos = 0;
  for(int64_t part=0; part<nparticles; ++part){
    
    while(!in_radius(xposition,yposition,x_c,y_c,radius)){
      xposition+= ldistance; 
      if(xposition > maxxposition){
        if(yposition > maxyposition){
          break;
        }
        xposition=x_topproc;
        yposition+=ldistance;
      }
    }

    if(xposition > maxxposition){
      if(yposition > maxyposition){
          break;
      }
    }

    tparticles++;
    x[part] = xposition;
    y[part] = yposition;
         
    xposition+=ldistance;
    if(xposition > maxxposition){
      if(yposition > maxyposition){
        break;
      }
      xposition=x_topproc;
      yposition+=ldistance;
    }


    P[part] = pressure_1;
    rho[part] = rho_1; 
    //u[part] = u_1;
    m[part] = density(x[part],y[part],x_c,y_c);
    u[part] = u_1;///m[part];
    h[part] = smoothing_length;
    id[part] = posid++; 
  }

  std::cout<<"Real number of particles: "<<tparticles<<std::endl;

  char filename[128];
  //sprintf(filename,"%s_%d.h5part",fileprefix,nparticles);
  sprintf(filename,"%s.h5part",fileprefix);
  // Remove the previous file 
  remove(filename); 


  Flecsi_Sim_IO::HDF5ParticleIO testDataSet; 
  testDataSet.createDataset(filename,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",tparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",dimension);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  //const char * simName = "sodtube_1D";
  //testDataSet.writeDatasetAttributeArray("name","string",simName);
  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",tparticles,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",tparticles,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",tparticles,z);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",tparticles,vx);
  _d2.createVariable("vy",Flecsi_Sim_IO::point,"double",tparticles,vy);
  _d3.createVariable("vz",Flecsi_Sim_IO::point,"double",tparticles,vz);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("ax",Flecsi_Sim_IO::point,"double",tparticles,ax);
  _d2.createVariable("ay",Flecsi_Sim_IO::point,"double",tparticles,ay);
  _d3.createVariable("az",Flecsi_Sim_IO::point,"double",tparticles,az);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();


  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",tparticles,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",tparticles,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",tparticles,u);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",tparticles,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",tparticles,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",tparticles,id);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  testDataSet.closeFile(); 
  
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] ax;
  delete[] ay;
  delete[] az;
  delete[] h;
  delete[] rho;
  delete[] u;
  delete[] P;
  delete[] m;
  delete[] id;
  delete[] dt;
 
  MPI_Finalize();
  return 0;
}
