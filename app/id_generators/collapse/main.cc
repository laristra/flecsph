/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>
#include <H5hut.h>

#include "user.h"
#include "collapse.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

/*
Cold Dust Cloud Collapse test
*/


//
// help message
//
void print_usage() {
  std::cout
      << "Initial data generator for the " << gdimension << "D Dust Cloud Collapse test"
      << std::endl << "Usage: ./collapse_generator <parameter-file.par>"
      << std::endl;
}

//
// derived parameters
//
static double timestep = 1.0;         // Recommended timestep
static double total_mass = 1.;        // total mass of the fluid
static double mass_particle = 1.;     // mass of an individual particle
static point_t bbox_max, bbox_min;    // bounding box of the domain
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  particle_lattice::select();

  // The value for constant timestep
  timestep = initial_dt;

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
  if (gdimension < 3) {
    clog(error) << "This test must be run in 3D" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
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
    clog(error) << "ERROR: parameter file not specified!" << std::endl;
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

  // Assign density, pressure and specific internal energy to particles
  for(int64_t part=0; part<nparticles; ++part){
    m[part] = mass_particle;
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    u[part] = uint_initial;
    h[part] = sph_smoothing_length;
    id[part] = posid++;
  }

  clog(info) << "Number of particles: " << nparticles << std::endl;

  // remove the previous file
  remove(initial_data_file.c_str());

  h5_file_t * dataFile = H5OpenFile(initial_data_file.c_str()
      ,H5_O_WRONLY, MPI_COMM_WORLD);
    
  int use_fixed_timestep = 1; 
  // add the global attributes
  
  H5WriteFileAttribInt64(dataFile,"nparticles",&nparticles,1);
  H5WriteFileAttribFloat64(dataFile,"timestep",&timestep,1);
  int dim = gdimension; 
  H5WriteFileAttribInt32(dataFile,"dimension",&dim,1);
  H5WriteFileAttribInt32(dataFile,"use_fixed_timestep",&use_fixed_timestep,1);

  H5SetStep(dataFile,0);
  H5PartSetNumParticles(dataFile,nparticles);
  H5PartWriteDataFloat64(dataFile,"x",x);
  H5PartWriteDataFloat64(dataFile,"y",y);
  H5PartWriteDataFloat64(dataFile,"z",z);
  H5PartWriteDataFloat64(dataFile,"vx",vx);
  H5PartWriteDataFloat64(dataFile,"vy",vy);
  H5PartWriteDataFloat64(dataFile,"vz",vz);
  H5PartWriteDataFloat64(dataFile,"ax",ax);
  H5PartWriteDataFloat64(dataFile,"ay",ay);
  H5PartWriteDataFloat64(dataFile,"az",az);
  H5PartWriteDataFloat64(dataFile,"h",h);
  H5PartWriteDataFloat64(dataFile,"rho",rho);
  H5PartWriteDataFloat64(dataFile,"u",u);
  H5PartWriteDataFloat64(dataFile,"P",P);
  H5PartWriteDataFloat64(dataFile,"m",m);
  H5PartWriteDataInt64(dataFile,"id",id);
 
  H5CloseFile(dataFile);

  delete[] x, y, z, vx, vy, vz, ax, ay, az, h, rho, u, P, m, id, dt; 

  MPI_Finalize();
  return 0;
}

