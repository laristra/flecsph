/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
#include <math.h>

#include "user.h"
#include "noh.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"
#include "io.h"
using namespace io;

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))


//
// help message
//
void print_usage() {
  clog_one(warn)
      << "Initial data generator for the " << gdimension << "D Noh collapse test"
      << std::endl << "Usage: ./noh_generator <parameter-file.par>" << std::endl;
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
  bbox_min = -sphere_radius;
  bbox_max =  sphere_radius;

  // particle separation
  SET_PARAM(sph_separation, (2.*sphere_radius/(lattice_nx-1)));

  // Count number of particles
  int64_t tparticles =
      particle_lattice::count(lattice_type,domain_type,
      bbox_min,bbox_max,sph_separation,0);
  SET_PARAM(nparticles, tparticles);

  // total mass
  if (gdimension == 2) {
    total_mass = rho_initial * M_PI*SQ(sphere_radius);
  }
  else if (gdimension == 3) {
    total_mass = rho_initial * 4./3.*M_PI*CU(sphere_radius);
  }
  else
    assert (false);

  // single particle mass
  assert (equal_mass);
  mass_particle = total_mass / nparticles;

  // set kernel
  kernels::select();

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

  clog_one(info) << "Sphere: r=" << sphere_radius << std::endl
                 << "origin: pos=["<<0<<";"<<0<<"]" << std::endl
                 << "Generating "  << nparticles << " particles in "
                 << gdimension << "D" << std::endl;

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

  // Assign density, pressure and specific internal energy to particles, etc.
  for(int64_t part=0; part<nparticles; ++part){
    m[part] = mass_particle;
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    ax[part] = 0.0;
    ay[part] = 0.0;
    az[part] = 0.0;
    u[part] = uint_initial;
    h[part] = sph_smoothing_length;
    id[part] = posid++;

    // Assign particle inward pointing velocity with absolute value 0.1
    double A = sqrt(SQ(x[part]) + SQ(y[part]) + SQ(z[part]));
    if (A <= 0.0) {
      vx[part] = 0.0;
      vy[part] = 0.0;
      if(gdimension > 2) vz[part] = 0.0;
    }
    else {
      vx[part] = -(x[part]) * 0.1 / A;
      vy[part] = -(y[part]) * 0.1 / A;
      if(gdimension > 2) vz[part] = -(z[part]) * 0.1 / A;
    }
  }

  // remove the previous file
  remove(initial_data_file.c_str());
  hid_t dataFile = H5P_openFile(initial_data_file.c_str(),H5F_ACC_RDWR);

  int use_fixed_timestep = 1;
  // add the global attributes
  H5P_writeAttribute(dataFile,"nparticles",&nparticles);
  H5P_writeAttribute(dataFile,"timestep",&timestep);
  int dim = gdimension;
  H5P_writeAttribute(dataFile,"dimension",&dim);
  H5P_writeAttribute(dataFile,"use_fixed_timestep",&use_fixed_timestep);

  H5P_setNumParticles(nparticles);
  H5P_setStep(dataFile,0);

  //H5PartSetNumParticles(dataFile,nparticles);
  H5P_writeDataset(dataFile,"x",x,nparticles);
  H5P_writeDataset(dataFile,"y",y,nparticles);
  H5P_writeDataset(dataFile,"z",z,nparticles);
  H5P_writeDataset(dataFile,"vx",vx,nparticles);
  H5P_writeDataset(dataFile,"vy",vy,nparticles);
  H5P_writeDataset(dataFile,"vz",vz,nparticles);
  H5P_writeDataset(dataFile,"h",h,nparticles);
  H5P_writeDataset(dataFile,"rho",rho,nparticles);
  H5P_writeDataset(dataFile,"u",u,nparticles);
  H5P_writeDataset(dataFile,"P",P,nparticles);
  H5P_writeDataset(dataFile,"m",m,nparticles);
  H5P_writeDataset(dataFile,"id",id,nparticles);

  H5P_closeFile(dataFile);
  delete[] x, y, z, vx, vy, vz, ax, ay, az, h, rho, u, P, m, id, dt;

  MPI_Finalize();
  return 0;
}
