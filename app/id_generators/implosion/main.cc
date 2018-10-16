/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>

#include "user.h"
#include "implosion.h"
#include "params.h"
#include "hdf5ParticleIO.h"
#include "lattice.h"
#include "kernels.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

/*
Implosion of a homogeneous disk. For details, see
I. Sagert, W.P. Even, T.T. Strother, PHYSICAL REVIEW E 95, 053206 (2017)
] D. García-Senz, A. Relaño, R. M. Cabezón, and E. Bravo, Mon. Not. R.
Astron. Soc. 392, 346 (2009).
*/


//
// help message
//
void print_usage() {
  clog(warn)
      << "Initial data generator for the " << gdimension << "D Implosion"
      << std::endl << "Usage: ./implosion_generator <parameter-file.par>"
      << std::endl;
}

//
// derived parameters
//
static double timestep = 1.0;         // Recommended timestep
static double inner_radius = 0.8;     // Radius of injection region
static double total_mass = 1.;        // total mass of the fluid
static double mass_particle = 1.;     // mass of an individual particle
static point_t bbox_max, bbox_min;    // bounding box of the domain
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  particle_lattice::select();

  // The value for constant timestep
  timestep = initial_dt;

  inner_radius = 0.8*sphere_radius;
  
  // Bounding box of the domain
  if (domain_type == 0) { // box
    bbox_min[0] = -box_length/2.;
    bbox_max[0] =  box_length/2.;
    bbox_min[1] = -box_width/2.;
    bbox_max[1] =  box_width/2.;
    if (gdimension == 3) {
      bbox_min[2] = -box_height/2.;
      bbox_max[2] =  box_height/2.;
    }
  }
  else if (domain_type == 1) { // sphere or circle
    bbox_min = -sphere_radius;
    bbox_max =  sphere_radius;
  }

  // particle separation
  if (domain_type == 0) {
    SET_PARAM(sph_separation, (box_length/(lattice_nx-1)));
  }
  else if (domain_type == 1) {
    SET_PARAM(sph_separation, (2.*sphere_radius/(lattice_nx-1)));
  }

  // Count number of particles
  int64_t tparticles =
      particle_lattice::count(lattice_type,domain_type,
      bbox_min,bbox_max,sph_separation,0);
  SET_PARAM(nparticles, tparticles);

  // total mass
  if (gdimension == 2) {
    if (domain_type == 0) // a box
      total_mass = rho_initial * box_length*box_width;
    else if (domain_type == 1) // a circle
      total_mass = rho_initial * M_PI*SQ(sphere_radius);
  }
  else if (gdimension == 3) {
    if (domain_type == 0) // a box
      total_mass = rho_initial * box_length*box_width*box_height;
    else if (domain_type == 1) // a sphere
      total_mass = rho_initial * 4./3.*M_PI*CU(sphere_radius);
  }
  else
    assert (false);

  // single particle mass
  assert (equal_mass);
  mass_particle = total_mass / nparticles;

  // set kernel
  kernels::select(sph_kernel);

  // smoothing length
  const double sph_h = sph_eta * kernels::kernel_width
                               * pow(mass_particle/rho_initial,1./gdimension);
  SET_PARAM(sph_smoothing_length, sph_h);

  // intial internal energy
  SET_PARAM(uint_initial, (pressure_initial/(rho_initial*(poly_gamma-1.0))));


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
    clog_one(error) << "ERROR: parameter file not specified!" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // Initialize the arrays to be filled later
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

  // Generate the lattice
  assert(nparticles ==
      particle_lattice::generate(lattice_type,domain_type,
      bbox_min,bbox_max,sph_separation,0, x, y, z));

  // Particle id number
  int64_t posid = 0;

  // Number of particles in the blast zone
  int64_t particles_blast = 0;

  // Total mass of particles in the blast zone
  double mass_blast = 0.;

  double particle_radius = 0.;
  // Count the number of particles and mass in the blast zone
  // The blast is centered at the origin ({0,0} or {0,0,0})
  for(int64_t part=0; part<nparticles; ++part)
    particle_radius = sqrt(SQ(x[part]) + SQ(y[part]) + SQ(z[part]));
    if (particle_radius >= inner_radius) {
       particles_blast++;
       mass_blast += mass_particle;
    }

  // Assign density, pressure, etc. to particles
  for(int64_t part=0; part<nparticles; ++part){
    m[part] = mass_particle;
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    u[part] = uint_initial;
    h[part] = sph_smoothing_length;
    id[part] = posid++;

    particle_radius = sqrt(SQ(x[part]) + SQ(y[part]) + SQ(z[part]));
    if (particle_radius >= inner_radius) {
      double m = 1.0;
      double y = 1.0;
      u[part] = (inner_radius - particle_radius)*m + y;
      P[part] = (poly_gamma -1.0) * rho[part] * u[part];
    }
  }

  clog(info) << "Number of particles: " << nparticles << std::endl;
  clog(info) << "Total number of seeded blast particles: " << particles_blast << std::endl;
//  clog(info) << "Total blast energy (E_blast = u_blast * total mass): "
//                 << sedov_blast_energy * mass_blast << std::endl;
  // remove the previous file
  remove(initial_data_file.c_str());

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(initial_data_file.c_str(),MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",nparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",gdimension);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",nparticles,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",nparticles,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",nparticles,z);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",nparticles,vx);
  _d2.createVariable("vy",Flecsi_Sim_IO::point,"double",nparticles,vy);
  _d3.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticles,vz);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("ax",Flecsi_Sim_IO::point,"double",nparticles,ax);
  _d2.createVariable("ay",Flecsi_Sim_IO::point,"double",nparticles,ay);
  _d3.createVariable("az",Flecsi_Sim_IO::point,"double",nparticles,az);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();


  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",nparticles,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",nparticles,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",nparticles,u);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",nparticles,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",nparticles,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",nparticles,id);

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