/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>

#include "params.h"
#include "hdf5ParticleIO.h"
#include "kernels.h"

//
// help message
//
void print_usage() {
  using namespace std;
  clog_one(warn) << "Initial data generator for Sod shocktube test in 1D" << endl
         << "Usage: ./sodtube_generator <parameter-file.par>" << endl;
}

//
// derived parameters
//
static int64_t nparticlesproc;        // number of particles per proc
static double rho_1, rho_2;           // densities
static double vx_1, vx_2;             // velocities
static double pressure_1, pressure_2; // pressures
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params(int rank, int size) {
  using namespace std;
  using namespace param;

  // number of particles per core
  nparticlesproc = nparticles/size + 1;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  // particle spacing and smoothing length
  SET_PARAM(sph_separation, 1.0/(double)nparticles);
  SET_PARAM(sph_smoothing_length, (sph_separation*20)); // TODO: use sph_eta instead

  // test selector
  switch (sodtest_num) {
    case (1):
      // -- left side      | right side -- //
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
  set_derived_params(rank,size);

  // screen output
  clog_one(info) << "Sod test #" << sodtest_num << " in 1D:" << endl
         << " - number of particles: " << nparticles << endl
         << " - particles per core:  " << nparticlesproc << endl
         << " - generated initial data file: "
         << initial_data_file << endl;

  // allocate arrays

  // Position
  double* x = new double[nparticlesproc]();
  double* y = new double[nparticlesproc]();
  double* z = new double[nparticlesproc]();
  // Velocity
  double* vx = new double[nparticlesproc]();
  double* vy = new double[nparticlesproc]();
  double* vz = new double[nparticlesproc]();
  // Acceleration
  double* ax = new double[nparticlesproc]();
  double* ay = new double[nparticlesproc]();
  double* az = new double[nparticlesproc]();
  // Smoothing length
  double* h = new double[nparticlesproc]();
  // Density
  double* rho = new double[nparticlesproc]();
  // Internal Energy
  double* u = new double[nparticlesproc]();
  // Pressure
  double* P = new double[nparticlesproc]();
  // Mass
  double* m = new double[nparticlesproc]();
  // Id
  int64_t* id = new int64_t[nparticlesproc]();
  // Timestep
  double* dt = new double[nparticlesproc]();

  // Generate data
  // Find middle to switch m, u and rho
  double middle = nparticles*sph_separation/2.;
  // Find my first particle position
  double lposition = sph_separation*nparticlesproc*rank;
  // Id of my first particle
  int64_t posid = nparticlesproc*rank;

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*max(pressure_1/rho_1,pressure_2/rho_2));

  // The value for constant timestep
  double timestep = 0.5*sph_separation/cs;


  for(int64_t part=0; part<nparticlesproc; ++part){
    id[part] = posid++;
    x[part] = lposition;

    if(x[part] > middle){
      P[part] = pressure_2;
      rho[part] = rho_2;
      vx[part] = vx_2;
    }else{
      P[part] = pressure_1;
      rho[part] = rho_1;
      vx[part] = vx_1;
    }

    // compute internal energy using gamma-law eos
    u[part] = P[part]/(poly_gamma-1.)/rho[part];

    // particle masses and smoothing length
    m[part] = rho[part]*sph_smoothing_length/10.;
    h[part] = sph_smoothing_length;

    // P,Y,Z,VY,VZ,AX,AY,AZ stay 0
    // Move to the next particle
    lposition += sph_separation;

  } // for part=0..nparticles

  // delete the output file if exists
  remove(initial_data_file.c_str());

  // Header data
  // the number of particles = nparticles
  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(initial_data_file,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",nparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",1);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  //testDataSet.writeDatasetAttributeArray("name","string",simName);
  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",nparticlesproc,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",nparticlesproc,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",nparticlesproc,z);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",nparticlesproc,vx);
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


  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",nparticlesproc,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",nparticlesproc,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",nparticlesproc,u);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",nparticlesproc,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",nparticlesproc,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",nparticlesproc,id);

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
