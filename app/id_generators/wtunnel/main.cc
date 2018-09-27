/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>

#include "user.h"
#include "sodtube.h"
#include "params.h"
#include "hdf5ParticleIO.h"
#include "lattice.h"

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
 * aerodynamical properties.
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
  int64_t parts_mid= 0;
  int64_t parts_lr = 0;
  double mass = 0;
  bool equal_separation = !equal_mass;
  tparticles =  particle_lattice::count(lattice_type,2,cbox_min,cbox_max,
                                        sph_separation,0);

  double lr_sph_sep = 0.;
  double temp_part = 0;
  double temp_part_new = 0;

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

  tparticles =  particle_lattice::generate(lattice_type,2,cbox_min,cbox_max,
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

  // Header data
  // the number of particles = nparticles
  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(initial_data_file,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",tparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",gdimension);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  //testDataSet.writeDatasetAttributeArray("name","string",simName);
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
  //_d2.createVariable("vy",Flecsi_Sim_IO::point,"double",nparticlesproc,vy);
  //_d3.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticlesproc,vz);

  testDataSet.vars.push_back(_d1);
  //testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  //_d1.createVariable("ax",Flecsi_Sim_IO::point,"double",nparticlesproc,ax);
  //_d2.createVariable("ay",Flecsi_Sim_IO::point,"double",nparticlesproc,ay);
  //_d3.createVariable("az",Flecsi_Sim_IO::point,"double",nparticlesproc,az);

  //testDataSet.vars.push_back(_d1);
  //testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  //testDataSet.writeVariables();


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
