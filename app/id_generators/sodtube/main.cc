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

//
// help message
//
void print_usage() {
  using namespace std;
  clog_one(warn) << "Initial data generator for Sod shocktube test in" <<
  gdimension << "D" << endl << "Usage: ./sodtube_generator <parameter-file.par>"
  << endl;
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

  // particle spacing and smoothing length
  SET_PARAM(sph_separation, (box_length/(double)(lattice_nx - 1)));
  if(gdimension==1){
    SET_PARAM(sph_smoothing_length, (sph_separation*25.)); // TODO: use sph_eta instead
  } else if(gdimension==2){
    SET_PARAM(sph_smoothing_length, (sph_separation*4.)); // TODO: ???
  } else if(gdimension==3){
    SET_PARAM(sph_smoothing_length, (sph_separation*3.)); // TODO: ???
  }

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
      clog_one(error) << "ERROR: invalid test (" << sodtest_num << ")." << endl;
      MPI_Finalize();
      exit(-1);
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
  clog_one(info) << "Sod test #" << sodtest_num << " in " << gdimension
         << "D:" << endl << " - number of particles: " << nparticles
         << endl << " - particles per core:  " << nparticlesproc << endl
         << " - generated initial data file: " << initial_data_file << endl;

  // allocate arrays
  int64_t tparticles = 0;
  int64_t parts_mid= 0;
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
  double cs = sqrt(poly_gamma*max(pressure_1/rho_1,pressure_2/rho_2));

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
      h[part] = sph_smoothing_length;


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
      h[part] = sph_smoothing_length;

    } // for part=0..nparticles
  }
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
