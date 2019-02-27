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
#include "io.h"
using namespace io;

/*
 * Wind tunnel test
 * ----------------
 * In the wind tunnel problem, a flow is confined to a square well potential in
 * y and z directions. It is initialized in a section of the tunnel upstream with
 * the following parameters:
 *  - box_length:             length of the section which contains intial flow;
 *  - box_width, box_height:  the size of the yz-well in y- and z-directions;
 *  - flow_velocity:          initial velocity;
 *  - rho_initial, pressure_initial
 * Different obstacles (e.g. airfoil) can be placed in the tunnel to study their
 * aerodynamic properties.
 *
 */
//
// help message
//
void print_usage() {
  using namespace std;
  clog_one(warn) << "Initial data generator for the wind tunnel test in"
                 << gdimension << "D" << endl << "Usage: ./wtunnel_"
                 << gdimension << "d_generator <parameter-file.par>"
                 << endl;
}

//
// derived parameters
//
static int64_t nparticlesproc;        // number of particles per proc
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

// geometric extents of the flow (box-shaped)
static point_t cbox_min, cbox_max;

void set_derived_params() {
  using namespace std;
  using namespace param;
  assert (gdimension > 1);

  // compute the total number of particles
  int64_t npd = lattice_nx;

  // x-dimension
  cbox_min[0] = 0.5*box_width;
  cbox_max[0] = cbox_min[0] + box_length;

  // y-dimension
  cbox_max[1] = box_width/2.0;
  cbox_min[1] =-box_width/2.0;
  npd *= (int64_t)((double)lattice_nx*box_width/box_length);

  // 3D case
  if (gdimension>2) {
     cbox_max[2] = box_height/2.0;
     cbox_min[2] =-box_height/2.0;
     npd *= (int64_t)((double)lattice_nx*box_height/box_length);
  }
  SET_PARAM(nparticles, npd);

  // particle spacing and smoothing length
  SET_PARAM(sph_separation, (box_length/(double)(lattice_nx - 1)));
  if(gdimension==2){
    SET_PARAM(sph_smoothing_length, (sph_separation*4.)); // TODO: ???
  } else if(gdimension==3){
    SET_PARAM(sph_smoothing_length, (sph_separation*3.)); // TODO: ???
  }

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace std;
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

  // screen output
  clog_one(info) << "Wind tunnel problem in " << gdimension
         << "D:" << endl << " - number of particles: " << nparticles
         << endl << " - particles per core:  " << nparticlesproc << endl
         << " - generated initial data file: " << initial_data_file << endl;

  // allocate arrays
  int64_t tparticles = 0;
  tparticles =  particle_lattice::count(lattice_type,0,cbox_min,cbox_max,
                                        sph_separation,0);

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

  tparticles =  particle_lattice::generate(lattice_type,0,cbox_min,cbox_max,
                                           sph_separation,0,x,y,z);

  // particle id number
  int64_t posid = 0;

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*pressure_initial/rho_initial);

  // The value for constant timestep
  double timestep = 0.5*sph_separation/cs;

  for(int64_t part=0; part<tparticles; ++part){
    id[part] = posid++;
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    vx[part] = -flow_velocity;
    m[part] = rho_initial/(double)tparticles;

    // compute internal energy using gamma-law eos
    u[part] = pressure_initial/(poly_gamma-1.)/rho_initial;

    // particle smoothing length
    h[part] = sph_smoothing_length;

  } // for part=0..nparticles

  clog_one(info) << "Actual number of particles: " << tparticles << std::endl;
  // delete the output file if exists
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
