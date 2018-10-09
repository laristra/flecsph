/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>

#include "user.h"
#include "collapse.h"
#include "params.h"
#include "hdf5ParticleIO.h"
#include "lattice.h"

/*
Homologous dust collapse for self-gravitating problems with spherical
flow geometry and negligible pressure. Test was described in Colage 
and White (1966) and used in e.g. Moenchmeyer and Mueller (1989) and 
as DustCollapse problem in the FLASH code
*/


//
// help message
//
void print_usage() {
  std::cout
      << "Initial data generator for the " << gdimension << "D homologous collapse"
      << std::endl << "Usage: ./collapse_generator <parameter-file.par>"
      << std::endl;
}

//
// derived parameters
//
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  // total number of particles
  if (gdimension < 3) {
    clog(error) << "This test must be run in 3D" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
  } 
  else {
    SET_PARAM(nparticles, lattice_nx*lattice_nx*lattice_nx);
  }

  SET_PARAM(uint_initial, (pressure_initial/(rho_initial*(poly_gamma-1.0))));
  SET_PARAM(sph_smoothing_length, (2.5*sph_separation));

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}

int main(int argc, char * argv[]){
  using namespace param;

  // launch MPI
  int rank, size;
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // check options list: exactly one option is allowed
  if (argc != 2) {
    clog(error) << "ERROR: parameter file not specified!" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  particle_lattice::select();

  // Header data
  // the number of particles = nparticles
  // The value for constant timestep
  double timestep = initial_dt;
  double radius   = sph_separation*(lattice_nx-1.)/2.;
  const point_t bbox_max =  radius;
  const point_t bbox_min = -radius;

  // Central coordinates: in most cases this should be centered at (0,0,0)
  double x_c = 0.0;
  double y_c = 0.0;
  double z_c = 0.0;

  // Count number of particles
  int64_t tparticles = 
      particle_lattice::count(lattice_type,domain_type,
      bbox_min,bbox_max,sph_separation,0); 

  // Initialize the arrays to be filled later
  // Position
  double* x = new double[tparticles]();
  double* y = new double[tparticles]();
  double* z = new double[tparticles]();
  // Velocity
  double* vx = new double[tparticles]();
  double* vy = new double[tparticles]();
  double* vz = new double[tparticles]();
  // Acceleration
  double* ax = new double[tparticles]();
  double* ay = new double[tparticles]();
  double* az = new double[tparticles]();
  // Smoothing length
  double* h = new double[tparticles]();
  // Density
  double* rho = new double[tparticles]();
  // Internal Energy
  double* u = new double[tparticles]();
  // Pressure
  double* P = new double[tparticles]();
  // Mass
  double* m = new double[tparticles]();
  // Id
  int64_t* id = new int64_t[tparticles]();
  // Timestep
  double* dt = new double[tparticles]();

  // Generate the lattice
  assert(tparticles ==
      particle_lattice::generate(lattice_type,domain_type,
      bbox_min,bbox_max,sph_separation,0, x, y, z));

  // particle id number
  int64_t posid = 0;

  // Number of particles in the blast zone
  int64_t particles_blast = 0;

  double mass;
  // Total mass of particles in the blast zone
  if(gdimension<3){
    clog(error) << "This test must be run in 3D" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
  } else{
    mass = rho_initial * 4./3.*M_PI*pow(radius,3.)/tparticles;
  }
  // Assign density, pressure, mass and specific internal energy to particles
  for(int64_t part=0; part<tparticles; ++part){
    m[part] = mass;
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    u[part] = uint_initial;
    h[part] = sph_smoothing_length;
    id[part] = posid++;
  }

  clog(info) << "Real number of particles: " << tparticles << std::endl;

  // remove the previous file
  remove(initial_data_file.c_str());

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(initial_data_file.c_str(),MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",tparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",gdimension);
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
