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
#include "implosion.h"
#include "params.h"
#include "lattice.h"

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
  clog_one(warn)
      << "Initial data generator for the " << gdimension << "D Implosion"
      << std::endl << "Usage: ./implosion_generator <parameter-file.par>"
      << std::endl;
}

//
// derived parameters
//
static double r_blast;      // Radius of injection region
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  // total number of particles
  if (gdimension == 2) {
    SET_PARAM(nparticles, lattice_nx*lattice_nx);
  } 
  else if (gdimension == 3) {
    SET_PARAM(nparticles, lattice_nx*lattice_nx*lattice_nx);
  }
  else {
    clog_one(error) << "This test only works for 2D and 3D" << std::endl;    
    print_usage();
    MPI_Finalize();
    exit(0);
 }


  SET_PARAM(uint_initial, (pressure_initial/(rho_initial*(poly_gamma-1.0))));
  SET_PARAM(sph_smoothing_length, (5.*sph_separation));

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

  particle_lattice::select();

  // Header data
  // the number of particles = nparticles
  // The value for constant timestep
  double timestep = initial_dt;
  double radius   = sph_separation*(lattice_nx-1.)/2.;
  double inner_radius = 0.8*radius;
  const point_t bbox_max = radius;
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

  // particle mass, determined via the density and size of the disk/sphere
  double mass;
  if(gdimension==2){
    mass = rho_initial * M_PI*pow(radius,2.)/tparticles;
  } else{
    mass = rho_initial * 4./3.*M_PI*pow(radius,3.)/tparticles;
  }

  // Assign density, pressure, etc. to particles
  for(int64_t part=0; part<tparticles; ++part){
    m[part] = mass;
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    u[part] = uint_initial;
    h[part] = sph_smoothing_length;
    id[part] = posid++;

    double particle_radius = sqrt((x[part]-x_c)*(x[part]-x_c) + (y[part]-y_c)*(y[part]-y_c)
                             + (z[part]-z_c)*(z[part]-z_c));
    if (particle_radius >= inner_radius) {
      double m = 1.0;
      double y = 1.0;
      u[part] = (inner_radius - particle_radius)*m + y;
      P[part] = (poly_gamma -1) * rho[part] * u[part];
    }
  }

  clog_one(info) << "Real number of particles: " << tparticles << std::endl;

  // remove the previous file
  remove(initial_data_file.c_str());

  h5_file_t * dataFile = H5OpenFile(filename,H5_O_WRONLY, MPI_COMM_WORLD);
    
  int use_fixed_timestep = 1; 
  // add the global attributes
  H5WriteFileAttribInt64(dataFile,"nparticles",&nparticles,1);
  H5WriteFileAttribFloat64(dataFile,"timestep",&timestep,1);
  H5WriteFileAttribInt32(dataFile,"dimension",&dimension,1);
  H5WriteFileAttribInt32(dataFile,"use_fixed_timestep",&use_fixed_timestep,1);

  H5SetStep(dataFile,0);
  H5PartSetNumParticles(dataFile,nparticlesproc);
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
