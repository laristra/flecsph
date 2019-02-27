/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>

#include "user.h"
#include "sodtube.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"
#include "io.h"
using namespace io;

//
// help message
//
void print_usage() {
  clog_one(warn)
      << "Initial data generator for Sod shocktube test in"
      << gdimension << "D" << std::endl
      << "Usage: ./sodtube_generator <parameter-file.par>" << std::endl;
}

//
// derived parameters
//
static int64_t nparticlesproc;        // number of particles per proc
static double rho_1, rho_2;           // densities
static double vx_1, vx_2;             // velocities
static double pressure_1, pressure_2; // pressures
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

// geometric extents of the three regions: left, right and central
static point_t cbox_min, cbox_max;
static point_t rbox_min, rbox_max;
static point_t lbox_min, lbox_max;

void set_derived_params() {
  using namespace std;
  using namespace param;

  // compute the total number of particles
  int64_t npd = lattice_nx;

  // 1D setup
  cbox_max[0] =  box_length/6.;
  cbox_min[0] = -cbox_max[0];
  lbox_min[0] = -box_length/2.;
  lbox_max[0] =  cbox_min[0];
  rbox_min[0] =  cbox_max[0];
  rbox_max[0] = -lbox_min[0];

  // 2D case
  if (gdimension>1) {
     cbox_max[1] = lbox_max[1] = rbox_max[1] = box_width/2.0;
     cbox_min[1] = lbox_min[1] = rbox_min[1] =-box_width/2.0;
     npd *= (int64_t)((double)lattice_nx*box_width/box_length);
  }

  // 3D case
  if (gdimension>2) {
     cbox_max[2] = lbox_max[2] = rbox_max[2] = box_height/2.0;
     cbox_min[2] = lbox_min[2] = rbox_min[2] =-box_height/2.0;
     npd *= (int64_t)((double)lattice_nx*box_height/box_length);
  }
  SET_PARAM(nparticles, npd);

  // test selector
  switch (sodtest_num) {
    case (1):
      // -- middle         | left and right side -- //
      rho_1      = 1.0;      rho_2      = 0.125;
      pressure_1 = 1.0;      pressure_2 = 0.1;
      vx_1       = 0.0;      vx_2       = 0.0;
      break;

    case (2):
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 0.4;      pressure_2 = 0.4;
      vx_1       =-2.0;      vx_2       = 2.0;
      break;

    case (3):
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 1000.;    pressure_2 = 0.01;
      vx_1       = 0.0;      vx_2       = 0.0;
      break;

    case (4):
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 0.01;     pressure_2 = 100.;
      vx_1       = 0.0;      vx_2       = 0.0;
      break;

    case (5):
      rho_1      = 5.99924;  rho_2      = 5.99242;
      pressure_1 = 460.894;  pressure_2 = 46.0950;
      vx_1       = 19.5975;  vx_2       =-6.19633;
      break;

    case (6): // 1D equivalent to the Noh problem
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 1.e-6;    pressure_2 = 1.e-6;
      vx_1       = 1.0;      vx_2       =-1.0;
      break;

    default:
      clog_one(error) << "ERROR: invalid test (" << sodtest_num << ")." << endl;
      MPI_Finalize();
      exit(-1);

  }

  // particle spacing
  SET_PARAM(sph_separation, (box_length/(double)(lattice_nx - 1)));

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace param;

  // launch MPI
  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // check options list: exactly one option is allowed
  if (argc != 2) {
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();
  particle_lattice::select();

  // set kernel
  kernels::select();

  // screen output
  std::cout << "Sod test #" << sodtest_num << " in " << gdimension
       << "D:" << std::endl <<
       " - generated initial data file: " << initial_data_file << std::endl;

  // allocate arrays
  int64_t tparticles = 0;
  int64_t parts_mid = 0;
  int64_t parts_lr = 0;
  double mass = 0;
  bool equal_separation = !equal_mass;
  if(equal_separation){
    tparticles =  particle_lattice::count(lattice_type,2,cbox_min,cbox_max,
                                          sph_separation,0);
    parts_mid = tparticles;
    tparticles += particle_lattice::count(lattice_type,2,rbox_min,rbox_max,
                                          sph_separation,tparticles-1);
    parts_lr = tparticles-parts_mid;
    tparticles += particle_lattice::count(lattice_type,2,lbox_min,lbox_max,
                                          sph_separation,tparticles-1);
  }

  double lr_sph_sep = 0.;
  double temp_part = 0;
  double temp_part_new = 0;
  if(equal_mass){
    tparticles = particle_lattice::count(lattice_type,2,cbox_min,cbox_max,
                                         sph_separation,0);
    if(gdimension==1){
      mass = rho_1*sph_separation;
      lr_sph_sep = mass/rho_2;
    } else if(gdimension==2){
      mass = rho_1*sph_separation*sph_separation;
      if (lattice_type == 1 or lattice_type == 2)
        mass *= sqrt(3.0)/2.0;
      lr_sph_sep = sph_separation * sqrt(rho_1/rho_2);
    } else{
      mass = rho_1*sph_separation*sph_separation*sph_separation;
      if (lattice_type == 1 or lattice_type == 2)
        mass *= 1.0/sqrt(2.0);
      lr_sph_sep = sph_separation * cbrt(rho_1/rho_2);
    }
    rbox_min += (lr_sph_sep - sph_separation) / 2.0; // adjust rbox
    tparticles += particle_lattice::count(lattice_type,2,rbox_min,rbox_max,
                                          lr_sph_sep,tparticles);
    tparticles += particle_lattice::count(lattice_type,2,lbox_min,lbox_max,
                                          lr_sph_sep,tparticles);
  }

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

  if(equal_separation){
    tparticles =  particle_lattice::generate(lattice_type,2,cbox_min,cbox_max,
                                             sph_separation,0,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,rbox_min,rbox_max,
                                             sph_separation,tparticles,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,lbox_min,lbox_max,
                                             sph_separation,tparticles,x,y,z);
  }
  if(equal_mass){
    tparticles =  particle_lattice::generate(lattice_type,2,cbox_min,cbox_max,
                                             sph_separation,0,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,rbox_min,rbox_max,
                                             lr_sph_sep,tparticles,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,lbox_min,lbox_max,
                                             lr_sph_sep,tparticles,x,y,z);
  }

  // particle id number
  int64_t posid = 0;

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*std::max(pressure_1/rho_1,pressure_2/rho_2));

  // The value for constant timestep
  double timestep = 0.5*sph_separation/cs;

  if(equal_separation){
    for(int64_t part=0; part<tparticles; ++part){
      id[part] = posid++;
      if (particle_lattice::in_domain_1d(x[part],
          cbox_min[0], cbox_max[0], domain_type)) {
        P[part] = pressure_1;
        rho[part] = rho_1;
        vx[part] = vx_1;
        m[part] = rho[part]/(double)parts_mid;
      }
      else {
        P[part] = pressure_2;
        rho[part] = rho_2;
        vx[part] = vx_2;
        m[part] = rho[part]/(double)parts_lr;
      }

      // compute internal energy using gamma-law eos
      u[part] = P[part]/(poly_gamma-1.)/rho[part];

      // particle smoothing length
      h[part] = sph_eta * kernels::kernel_width
                        * pow(m[part]/rho[part],1./gdimension);
    } // for part=0..nparticles
  }
  if(equal_mass){
    for(int64_t part=0; part<tparticles; ++part){
      id[part] = posid++;
      if (particle_lattice::in_domain_1d(x[part],
          cbox_min[0], cbox_max[0], domain_type)) {
        P[part] = pressure_1;
        rho[part] = rho_1;
        vx[part] = vx_1;
      }
      else {
        P[part] = pressure_2;
        rho[part] = rho_2;
        vx[part] = vx_2;
      }

      // compute internal energy using gamma-law eos
      u[part] = P[part]/(poly_gamma-1.)/rho[part];

      // particle masses and smoothing length
      m[part] = mass;
      h[part] = sph_eta * kernels::kernel_width
                        * pow(m[part]/rho[part],1./gdimension);

    } // for part=0..nparticles
  }
  std::cout << "Actual number of particles: " << tparticles << std::endl
    << std::flush;
  // delete the output file if exists
  remove(initial_data_file.c_str());
  hid_t dataFile = H5P_openFile(initial_data_file.c_str(),H5F_ACC_RDWR);

  int use_fixed_timestep = 1;
  // add the global attributes
  H5P_writeAttribute(dataFile,"nparticles",&tparticles);
  H5P_writeAttribute(dataFile,"timestep",&timestep);
  int dim = gdimension;
  H5P_writeAttribute(dataFile,"dimension",&dim);
  H5P_writeAttribute(dataFile,"use_fixed_timestep",&use_fixed_timestep);

  H5P_setNumParticles(tparticles);
  H5P_setStep(dataFile,0);

  //H5PartSetNumParticles(dataFile,nparticles);
  H5P_writeDataset(dataFile,"x",x,tparticles);
  H5P_writeDataset(dataFile,"y",y,tparticles);
  H5P_writeDataset(dataFile,"z",z,tparticles);
  H5P_writeDataset(dataFile,"vx",vx,tparticles);
  H5P_writeDataset(dataFile,"vy",vy,tparticles);
  H5P_writeDataset(dataFile,"h",h,tparticles);
  H5P_writeDataset(dataFile,"rho",rho,tparticles);
  H5P_writeDataset(dataFile,"u",u,tparticles);
  H5P_writeDataset(dataFile,"P",P,tparticles);
  H5P_writeDataset(dataFile,"m",m,tparticles);
  H5P_writeDataset(dataFile,"id",id,tparticles);

  H5P_closeFile(dataFile);

  delete[] x, y, z, vx, vy, vz, ax, ay, az, h, rho, u, P, m, id, dt;

  MPI_Finalize();
  return 0;
}
