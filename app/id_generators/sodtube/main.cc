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
#include "sodtube.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"

//
// help message
//
void print_usage() {
  clog(warn)
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

    default:
      clog(error) << "ERROR: invalid test (" << sodtest_num << ")." << endl;
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
  kernels::select(sph_kernel);

  // screen output
  std::cout << "Sod test #" << sodtest_num << " in " << gdimension
       << "D:" << std::endl << " - number of particles: " << nparticles
       << std::endl << " - particles per core:  " << nparticlesproc << std::endl
       << " - generated initial data file: " << initial_data_file << std::endl;

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
      mass = rho_1*(cbox_max[0]-cbox_min[0])/tparticles;
      if(lattice_type==0){
        temp_part = tparticles;
      } else if(lattice_type==1 || lattice_type==2){
        temp_part = tparticles/sqrt(2.);
      }
      temp_part_new = rho_2/rho_1*(temp_part);
      lr_sph_sep = 1./(temp_part_new-1.);
    } else if(gdimension==2){
      mass = rho_1*(cbox_max[0]-cbox_min[0])*(cbox_max[1]-cbox_min[1])/tparticles;
      if(lattice_type==0){
        temp_part = tparticles;
      } else if(lattice_type==1 || lattice_type==2){
        temp_part = tparticles/sqrt(2.);
      }
      temp_part_new = rho_2/rho_1*(temp_part);
      lr_sph_sep = 1./((int)sqrt(temp_part_new)-1.);
    } else{
      mass = rho_1*(cbox_max[0]-cbox_min[0])*(cbox_max[1]-cbox_min[1])*(cbox_max[2]-cbox_min[2])/tparticles;
      if(lattice_type==0){
        temp_part = tparticles;
      } else if(lattice_type==1 || lattice_type==2){
        temp_part = tparticles/sqrt(2.);
      }
      temp_part_new = rho_2/rho_1*(temp_part);
      lr_sph_sep = 1./((int)cbrt(temp_part_new)-1.);
    }
    tparticles += particle_lattice::count(lattice_type,2,rbox_min,rbox_max,
                                          lr_sph_sep,tparticles-1);
    tparticles += particle_lattice::count(lattice_type,2,lbox_min,lbox_max,
                                          lr_sph_sep,tparticles-1);
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
                                             sph_separation,tparticles-1,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,lbox_min,lbox_max,
                                             sph_separation,tparticles-1,x,y,z);
  }
  if(equal_mass){
    tparticles =  particle_lattice::generate(lattice_type,2,cbox_min,cbox_max,
                                             sph_separation,0,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,rbox_min,rbox_max,
                                             lr_sph_sep,tparticles-1,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,lbox_min,lbox_max,
                                             lr_sph_sep,tparticles-1,x,y,z);
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
      if(x[part] < cbox_min[0] || x[part] > cbox_max[0]){
        P[part] = pressure_2;
        rho[part] = rho_2;
        vx[part] = vx_2;
        m[part] = rho[part]/(double)parts_lr;
      } else{
        P[part] = pressure_1;
        rho[part] = rho_1;
        vx[part] = vx_1;
        m[part] = rho[part]/(double)parts_mid;
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
      if(x[part] < cbox_min[0] || x[part] > cbox_max[0]){
        P[part] = pressure_2;
        rho[part] = rho_2;
        vx[part] = vx_2;
      } else{
        P[part] = pressure_1;
        rho[part] = rho_1;
        vx[part] = vx_1;
      }

      // compute internal energy using gamma-law eos
      u[part] = P[part]/(poly_gamma-1.)/rho[part];

      // particle masses and smoothing length
      m[part] = mass;
      h[part] = sph_eta * kernels::kernel_width
                        * pow(m[part]/rho[part],1./gdimension);

    } // for part=0..nparticles
  }
  clog_one(info) << "Actual number of particles: " << tparticles << std::endl
    << std::flush;
  // delete the output file if exists
  remove(initial_data_file.c_str());

  h5_file_t * dataFile = H5OpenFile(initial_data_file.c_str(),
      H5_O_WRONLY | H5_VFD_MPIIO_IND, MPI_COMM_WORLD);
    
  int use_fixed_timestep = 1; 
  // add the global attributes
  H5WriteFileAttribInt64(dataFile,"nparticles",&tparticles,1);
  H5WriteFileAttribFloat64(dataFile,"timestep",&timestep,1);
  int dim = gdimension; 
  H5WriteFileAttribInt32(dataFile,"dimension",&dim,1);
  H5WriteFileAttribInt32(dataFile,"use_fixed_timestep",&use_fixed_timestep,1);

  H5SetStep(dataFile,0);
  H5PartSetNumParticles(dataFile,tparticles);
  H5PartWriteDataFloat64(dataFile,"x",x);
  H5PartWriteDataFloat64(dataFile,"y",y);
  H5PartWriteDataFloat64(dataFile,"z",z);
  H5PartWriteDataFloat64(dataFile,"vx",vx);
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
