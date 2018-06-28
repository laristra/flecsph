/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>

#include "hdf5ParticleIO.h"
#include "kernels.h"


namespace simulation_params {
  int64_t nparticles;        // global number of particles
  int64_t nparticlesproc;    // number of particles per proc
  double ldistance;          // particles spacing
  double localgamma;         // polytropic index
  double smoothing_length;   // constant smoothing length

  // test conditions for two sides of the domain
  int    sodtest_num;            // which Sod test to generate
  double rho_1, rho_2;           // densities
  double vx_1, vx_2;             // velocities
  double pressure_1, pressure_2; // pressures

  // output filename
  const char* fileprefix = "hdf5_sodtube";
  char output_filename[128];
}


//
// setup parameter defaults
//
void set_default_param(int rank, int size) {
  using namespace simulation_params;

  // number of particles
  nparticles = 1000;

  // equation of state parameters (one so far)
  localgamma = 1.4;//5./3.;

  // run Sod test 1 by default
  sodtest_num = 1;
}


//
// help message
//
void print_usage(int rank) {
  using namespace std;
  clog(warn) << "Initial data generator for Sod shocktube test in 1D" << endl
         << "Usage: ./sodtube_generator [OPTIONS]" << endl
         << " -h: this help" << endl
         << " -n <number of particles>" << endl
         << " -t <Sod test (integer from 1 to 5)>" << endl;
}


//
// option parser
//
void parse_command_line_options(int rank, int size, int argc, char* argv[]) {
  using namespace std;
  using namespace simulation_params;

  for (int i=1; i<argc; ++i)
    if (argv[i][0] == '-')
      switch(argv[i][1]) {
      case 'h':
        print_usage(rank);
        MPI_Finalize();
        exit(0);
        break;

      case 'n':
        nparticles = atoll(argv[++i]);
        nparticlesproc = nparticles/size;
        if(rank==size-1)
          nparticlesproc = nparticles - nparticlesproc*(size-1);
        break;

      case 't':
        sodtest_num = atoi(argv[++i]);
        assert(sodtest_num>=1 && sodtest_num<=6);
        break;

      default:
        clog(error) << "ERROR: unknown option '-" << argv[i][1] << "'" << endl;
        MPI_Finalize();
        exit(-1);

      } // switch argv[i][1]

}


//
// setup simulation parameters
//
void set_param(int rank, int size) {
  using namespace std;
  using namespace simulation_params;

  // number of particles per core
  nparticlesproc = nparticles/size + 1;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  // particle spacing and smoothing length
  ldistance = 1.0/(double)nparticles;
  smoothing_length = ldistance*10; // TODO: introduce \eta parameter

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
      clog(error) << "ERROR: invalid test (" << sodtest_num << ")." << endl;
      MPI_Finalize();
      exit(-1);
  }

  // output file
  sprintf(output_filename,"%s.h5part",fileprefix);

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace std;
  using namespace simulation_params;

  // launch MPI 
  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  // set simulation parameters
  set_default_param(rank,size);
  parse_command_line_options(rank,size,argc,argv);
  set_param(rank,size);

  // screen output
  clog(info) << "Sod test #" << sodtest_num << " in 1D:" << endl
         << " - number of particles: " << nparticles << endl
         << " - particles per core:  " << nparticlesproc << endl
         << " - output file: " << output_filename << endl;

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
  double middle = nparticles*ldistance/2.;
  // Find my first particle position
  double lposition = ldistance*nparticlesproc*rank;
  // Id of my first particle
  int64_t posid = nparticlesproc*rank;

  // max. value for the speed of sound
  double cs = sqrt(localgamma*max(pressure_1/rho_1,pressure_2/rho_2));

  // The value for constant timestep
  double timestep = 0.5*ldistance/cs;


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
    u[part] = P[part]/(localgamma-1.)/rho[part];

    // particle masses and smoothing length
    m[part] = rho[part]*smoothing_length/10.;
    h[part] = smoothing_length;

    // P,Y,Z,VY,VZ,AX,AY,AZ stay 0
    // Move to the next particle
    lposition += ldistance;

  } // for part=0..nparticles

  // delete the output file if exists
  remove(output_filename);

  // Header data
  // the number of particles = nparticles
  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(output_filename,MPI_COMM_WORLD);

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
