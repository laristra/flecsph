/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


/*
 * Noh Collapse test
 * -----------------
 * The tet is initialized as disk / sphere of particles with homogeneous density, 
 * vanishingly small pressure, and inward velocity. As particles move inwards, they
 * pile up at the center at forming a region with constant density that is dependent 
 * on gamma and the dimensionality of the problem. 
 * A standing shock front forms that moves outwards as more particles are piling up.
 *
 * For an analytic solution, see the code noh.f in /src/tools
 * For more information, see: 
 * Liska & Wendroff, 2003, SIAM J. Sci. Comput. 25, 995, Section 4.5
*/


#include <iostream>
#include <algorithm>
#include <cassert>

#include "hdf5ParticleIO.h"
#include "kernel.h"

    
const double ldistance = 0.001;                     // Distance between the particles 
const double localgamma = 5./3.;                    // Gamma for ideal gas eos
const double rho_in = 1;                            // Initial density 
const double P_in = 1.0e-6;                         // Initial pressure
const double smoothing_length = 2.0*ldistance;      
const char* fileprefix = "hdf5_noh";                // Output file prefix


// Is particle with x,y coordinates within a disk (center located at x0,y0) with radius r?
bool in_radius(
    double x,
    double y,
    double x0,
    double y0, 
    double r) {
    return (x-x0)*(x-x0)+(y-y0)*(y-y0) < r*r;
}


int main(int argc, char * argv[]){

  int64_t sparticles = 100;         // Default number of particles

  if (argc != 2) {
    printf("./noh_generator [square nParticles]\n");
    fprintf(stderr,"Generating default number of particles=%ld*%ld=%ld",
    sparticles,sparticles,sparticles*sparticles);
  }
  else sparticles = atoll(argv[1]);

  int rank, size; 
  int provided; 

  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE); 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  // Total number of initialized particles
  int64_t nparticles = sparticles*sparticles;

  // Radius of the disk given the number of particles 
  double radius = (ldistance*(sparticles - 1))/2.0; 

  // Center of disk
  double x_c = (sparticles - 1)*ldistance/2.0;
  double y_c = x_c;

  // Particle mass given the initial density and number of particles
  double m_in = rho_in * radius * radius * 4.0 / nparticles; 

  std::cout<<"Sphere: r="<<radius<<" pos=["<<x_c<<";"<<y_c<<"]"<<std::endl;

  //if(rank==0){
    printf("Generating %ld particles by %ldx%ld in sphere r=%.4f\n",
        nparticles,sparticles,sparticles,radius);
  //}


  // Start on  0 0
  double x_topproc = x_c - radius;
  double y_topproc = y_c - radius;

  double maxxposition = x_c + radius;
  double maxyposition = y_c + radius;

  // Particle positions   
  double* x = new double[nparticles]();
  double* y = new double[nparticles]();
  double* z = new double[nparticles]();

  // Particle velocities
  double* vx = new double[nparticles]();
  double* vy = new double[nparticles]();
  double* vz = new double[nparticles]();

  // Particle accelerations 
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

  // ID
  int64_t* id = new int64_t[nparticles]();

  // Timestep
  double* dt = new double[nparticles]();
  
  // ID of first particle 
  int64_t posid = 0;


  // Header data 
  // the number of particles = nparticles 
  // The value for constant timestep 
  double timestep = 0.001;
  int dimension = 2;
  
  std::cout << "top_X=" << x_topproc << " top_Y=" << y_topproc << " maxX=" << maxxposition << " maxY=" << maxyposition << std::endl;

  double xposition = x_topproc; 
  double yposition = y_topproc;
  int64_t tparticles = 0;


  // Loop over all particles and assign position, velocity etc. 
  for (int64_t part = 0; part < nparticles; ++part) {    

    // Checks if particle is in the disk and creates position coordinates 
    while (!in_radius(xposition, yposition, x_c, y_c, radius)) {
      xposition += ldistance; 
      if (xposition > maxxposition) {
        if (yposition > maxyposition) {break;}
        xposition  = x_topproc;
        yposition += ldistance;
      }
    }
    if (xposition > maxxposition) {
      if (yposition > maxyposition) {break;}
    }

    tparticles++;

    // Assign particle position
    x[part] = xposition;
    y[part] = yposition;
         
    xposition += ldistance;
    if (xposition > maxxposition) {
      if (yposition > maxyposition) {break;}
      xposition  = x_topproc;
      yposition += ldistance;
    }

    // Assign particle pressure
    P[part] = P_in;

    // Assign particle density 
    rho[part] = rho_in; 

    // Assign particle accelerations 
    ax[part] = 0.0;
    ay[part] = 0.0;
    az[part] = 0.0;

    // Assign particle mass
    m[part] = m_in;

    // Assign particle internal energy 
    u[part] = P[part]/(localgamma-1.)/rho[part];

    // Assign particle smoothing length
    h[part] = smoothing_length;

    // Assign particle ID
    id[part] = posid++;

    // Assign particle inward pointing velocity with absolute value 0.1
    double A = sqrt((x[part]-x_c)*(x[part]-x_c) + (y[part]-y_c)*(y[part]-y_c)); 
    if (A <= 0.0) {
      vx[part] = 0.0;
      vy[part] = 0.0;
    }
    else {
      vx[part] = -(x[part]-x_c) * 0.1 / A;
      vy[part] = -(y[part]-y_c) * 0.1 / A;  
    }
  }

  std::cout << "Real Number of Particles: " << tparticles << std::endl;

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
