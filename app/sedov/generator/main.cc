/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>

#include "hdf5ParticleIO.h"
#include "kernels.h"



/*
The Sedov test is set up with uniform density and vanishingly small pressure.
An explosion is initialized via a point-like deposition of energy E_blast in
the center of the simulation space. The resulting spherically symmetric shock 
waves moves outwards with a radial distance given by

r_shock = (E_blast * t^2 / (alpha * n_0))^(1/(2*d))

where alpha is a constant of the order one with its exact value given by the 
adiabatic index gamma. The peak density of the shock is given by: 

n_shock = n_0 ((gamma + 1)/(gamma -1)) 

The velocity of shocked matter has a radial dependence approximately ~ r/t. 
As the shock wave moves away from the center it leaves behind matter at 
vanishingly low density. With the pressure staying finite for r=0, the 
temperature grows and becomes infinitely large at the origin of the blast wave . 

Reference:
G.  Taylor, “The Formation of a Blast Wave by a Very Intense Explosion.  
I.  Theoretical  Discussion,” Royal Society of London Proceedings Series A
, vol. 201, pp. 159–174, Mar. 1950.
*/


const double ldistance = 0.001;  // Distance between the particles 
const double localgamma = 1.4;   // Set to fit value in default_physics.h
const double rho_in = 1;
const double pressure_in = 1.0e-7;
const double u_in = pressure_in/(rho_in*(localgamma - 1.0));
const double smoothing_length = 5.*ldistance;
const char* fileprefix = "hdf5_sedov";

const double u_blast = 1.0;             // Injected total blast energy
const double r_blast = 0.5*ldistance;   // Radius of injection region 


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


  int64_t sparticles = 101;
  int rank, size; 
  int provided; 
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE); 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  if (argc != 2) {
    clog_one(warn) << "WARNING: you have not specified sqrt of the number of particles!" 
           << std::endl << "Usage: ./sedov_generator [sParticles]" << std::endl
           << " - generates initial conditions with  sParticles^2 particles."
           << std::endl << "Generating with the default value: sparticles = " 
           << sparticles << std::endl;
  }else{
    sparticles = atoll(argv[1]);
    clog_one(info) << "Square root of the number of particles: sparticles = " 
           << sparticles << std::endl;
  }

  int64_t nparticles = sparticles*sparticles;
  clog_one(info) << "Generating " << nparticles << " particles" << std::endl;

  // Start on  0 0

  double maxxposition = (sparticles)*ldistance;
  double maxyposition = (sparticles)*ldistance;


  // Central coordinates 
  double x_c = maxxposition/2.0; 
  double y_c = maxyposition/2.0; 


  // Particle mass from number of particles and density 
  double mass = rho_in * maxxposition * maxyposition/nparticles;  

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

  // Number of particles in the blast zone
  int64_t particles_blast = 0;

  // Total mass of particles in the blast zone
  double mass_blast = 0;  

  double xposition = 0;
  int64_t tparticles = 0;
  double yposition = 0;

  
  for (int64_t part=0; part<nparticles; ++part) {

    tparticles++;
    x[part] = xposition;
    y[part] = yposition;
    m[part] = mass;

    // Count particles in the blast zone and sum their masses 
    if(sqrt((x[part]-x_c)*(x[part]-x_c)+(y[part]-y_c)*(y[part]-y_c)) < r_blast){
       particles_blast++;
       mass_blast += m[part];
    }

    xposition+= ldistance;
    if(xposition > maxxposition){
      xposition = 0.;
      yposition+=ldistance;
    }
  }
  
  // Assign density, pressure and specific internal energy to particles, 
  // including the particles in the blast zone
  for(int64_t part=0; part<nparticles; ++part){
    
    P[part] = pressure_in;
    rho[part] = rho_in; 
    u[part] = u_in;
    h[part] = smoothing_length;
    id[part] = posid++;
 
    if(sqrt((x[part]-x_c)*(x[part]-x_c)+(y[part]-y_c)*(y[part]-y_c)) < r_blast){
       u[part] = u_blast/particles_blast;
       P[part] = u[part]*rho[part]*(localgamma - 1.0);
    }
  }

  clog_one(info) << "Real number of particles: " << tparticles << std::endl;
  clog_one(info) << "Total blast energy (E_blast = u_blast * total mass): " << u_blast * mass_blast << std::endl;

  char filename[128];
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
