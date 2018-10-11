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
#include "sedov.h"
#include "params.h"
#include "lattice.h"

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


//
// help message
//
void print_usage() {
  std::cout
      << "Initial data generator for the " << gdimension << "D Sedov blast wave"
      << std::endl << "Usage: ./sedov_generator <parameter-file.par>"
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
  else {
    SET_PARAM(nparticles, lattice_nx*lattice_nx*lattice_nx);
  }

  SET_PARAM(uint_initial, (pressure_initial/(rho_initial*(poly_gamma-1.0))));
  SET_PARAM(sph_smoothing_length, (5.*sph_separation));
  r_blast = sedov_blast_radius * sph_separation;   // Radius of injection region

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

  // Number of particles in the blast zone
  int64_t particles_blast = 0;

  // Total mass of particles in the blast zone
  double mass_blast = 0;
  double mass;
  if(gdimension==2){
    mass = rho_initial * M_PI*pow(radius,2.)/tparticles;
  } else{
    mass = rho_initial * 4./3.*M_PI*pow(radius,3.)/tparticles;
  }
  // Assign density, pressure and specific internal energy to particles,
  // including the particles in the blast zone
  for(int64_t part=0; part<tparticles; ++part){
    // Particle mass from number of particles and density
    m[part] = mass;

    // Count particles in the blast zone and sum their masses
    if(sqrt((x[part]-x_c)*(x[part]-x_c)+(y[part]-y_c)*(y[part]-y_c)+(z[part]-z_c)*(z[part]-z_c)) < r_blast){
       particles_blast++;
       mass_blast += m[part];
    }
  }
  for(int64_t part=0; part<tparticles; ++part){
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    u[part] = uint_initial;
    h[part] = sph_smoothing_length;
    id[part] = posid++;

    if(sqrt((x[part]-x_c)*(x[part]-x_c)+(y[part]-y_c)*(y[part]-y_c)+(z[part]-z_c)*(z[part]-z_c)) < r_blast){
       u[part] = u[part]+sedov_blast_energy/particles_blast;
       P[part] = u[part]*rho[part]*(poly_gamma - 1.0);
    }
  }

  clog(info) << "Real number of particles: " << tparticles << std::endl;
  clog(info) << "Total number of seeded blast particles: " << particles_blast << std::endl;
  clog(info) << "Total blast energy (E_blast = u_blast * total mass): "
                 << sedov_blast_energy * mass_blast << std::endl;

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
