/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


/*
 * Noh Collapse test
 * -----------------
 * The test is initialized as a disk / sphere of particles with homogeneous 
 * density, vanishingly small pressure, and inward velocity v_r. As particles 
 * move inwards, they pile up at the center at forming a region with constant 
 * density that is dependent on gamma and the dimensionality of the problem.
 * A standing shock front forms that moves outwards as more particles are piling
 * up. Its radial distance grows as 
 *
 * r_shock (t) = 0.5 (gamma - 1) v_r t
 *
 * The density of infalling matter evolves as:
 *
 * n(r>=r_shock) = n_0 (1 + (v_r/r)t)^(d-1)
 *
 * The density of matter enclose by the shock is given by:
 *
 * n(r<r_rhock) = n_0 * ((gamma + 1)/(gamma - 1))^d
 *
 * where d gives the geometry of the system (1=planar, 2=cyl., 3=spherical)
 *
 * For more information, see:
 * Liska & Wendroff, 2003, SIAM J. Sci. Comput. 25, 995, Section 4.5
*/


#include <iostream>
#include <algorithm>
#include <cassert>

#include "params.h"
#include "hdf5ParticleIO.h"
#include "kernels.h"

//
// help message
//
void print_usage() {
  clog_one(warn)
      << "Initial data generator for the Noh collapse test" << std::endl
      << "Usage: ./noh_generator <parameter-file.par>"      << std::endl;
}


//
// derived parameters
//
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  // total number of particles
  SET_PARAM(nparticles, sqrt_nparticles*sqrt_nparticles);
  SET_PARAM(sph_smoothing_length, (2.0*sph_separation));

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}


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
  using namespace param;

  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // check options list: exactly one command-line argument is accepted
  if (argc != 2) {
    clog_one(error) << "ERROR: parameter file not specified!" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // Radius of the disk given the number of particles
  double radius = (sph_separation*(sqrt_nparticles - 1))/2.0;

  // Center of disk
  double x_c = (sqrt_nparticles - 1)*sph_separation/2.0;
  double y_c = x_c;


  clog_one(info) << "Sphere: r=" << radius << std::endl
                 << "origin: pos=["<<x_c<<";"<<y_c<<"]" << std::endl
                 << "Generating "  << nparticles
                 << " particles (="<< sqrt_nparticles<<"^2) " << std::endl;


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
  double timestep = initial_dt;  // TODO: replace with parameter initial_dt
  int dimension = 2;             // TODO: use global dimension parameter

  clog_one(info) << "top_X=" << x_topproc << " top_Y="  << y_topproc    << std::endl
                 << "maxX=" << maxxposition << " maxY=" << maxyposition << std::endl;

  double xposition = x_topproc;
  double yposition = y_topproc;
  int64_t tparticles = 0;
  double max_radius = 0;

  // Loop over all particles and assign position, velocity etc.
  for (int64_t part = 0; part < nparticles; ++part) {

    // Checks if particle is in the disk and creates position coordinates
    while (!in_radius(xposition, yposition, x_c, y_c, radius)) {
      xposition += sph_separation;
      if (xposition > maxxposition) {
        if (yposition > maxyposition) {break;}
        xposition  = x_topproc;
        yposition += sph_separation;
      }
    }
    if (xposition > maxxposition) {
      if (yposition > maxyposition) {break;}
    }

    tparticles++;

    // Assign particle position
    x[part] = xposition;
    y[part] = yposition;

    // Determine the maximum radial distance for a particle for 
    // later mass/density determination 
    if (sqrt((x[part]-x_c)*(x[part]-x_c) + (y[part]-y_c)*(y[part]-y_c)) 
    > max_radius) {
        max_radius = sqrt((x[part]-x_c)*(x[part]-x_c) 
                        + (y[part]-y_c)*(y[part]-y_c));
    }

    xposition += sph_separation;
    if (xposition > maxxposition) {
      if (yposition > maxyposition) {break;}
      xposition  = x_topproc;
      yposition += sph_separation;
    }

    // Assign particle pressure
    P[part] = pressure_initial;

    // Assign particle density
    rho[part] = rho_initial;

    // Assign particle accelerations
    ax[part] = 0.0;
    ay[part] = 0.0;
    az[part] = 0.0;

    // Assign particle mass
    //m[part] = m_in;

    // Assign particle internal energy
    u[part] = P[part]/(poly_gamma-1.)/rho[part];

    // Assign particle smoothing length
    h[part] = sph_smoothing_length;

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

  clog_one(info) << "Actual number of particles inside the sphere: "
                 << tparticles << std::endl;

  // Particle mass given the density, number of particles, and disk radius
  double m_in = rho_initial * max_radius * max_radius * M_PI / tparticles;

  for (int64_t part = 0; part < nparticles; ++part) m[part] = m_in;


  // remove the previous file
  remove(initial_data_file.c_str());

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(initial_data_file.c_str(),MPI_COMM_WORLD);

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
