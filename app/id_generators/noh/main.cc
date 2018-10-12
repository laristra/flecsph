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
#include <math.h>

#include "user.h"
#include "noh.h"
#include "params.h"
#include "hdf5ParticleIO.h"
#include "lattice.h"
#include "kernels.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

//
// help message
//
void print_usage() {
  clog(warn)
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

  clog(info) << "Sphere: r=" << sphere_radius << std::endl
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
