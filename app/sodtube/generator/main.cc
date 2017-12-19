/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>

#include "hdf5ParticleIO.h"
#include "physics/kernel.h"


namespace simulation_params {
  int64_t nparticles;        // global number of particles
  int64_t nparticlesproc;    // number of particles per proc
  double ldistance;          // particles spacing 
  double localgamma;         // polytropic index
  double smoothing_length;   // constant smoothing length

  // test conditions for two sides of the domain
  int    sodtest_num;            // which Sod test to generate
  double rho_1, rho_2;           // densities
  double pressure_1, pressure_2; // pressures
  double u_1, u_2;               // internal energies
  double m_1, m_2;               // particle masses

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
  localgamma = 5./3.;

  // run Sod test 1 by default
  sodtest_num = 1;
}


//
// help message
//
void print_usage(int rank) {
  using namespace std;
  if (rank == 0)
    cout << "Initial data generator for Sod shocktube test in 1D" << endl
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
        cerr << "ERROR: unknown option '-" << argv[i][1] << "'" << endl;
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
      // -- left side   || right side -- //  
      rho_1      = 1.0;    rho_2      = 0.125;
      pressure_1 = 1.0;    pressure_2 = 0.1;
      u_1        = 2.5;    u_2        = 2.0;
      m_1        = 1.e-4;  m_2        = 1.e-5;
      break;

    default:
      cerr << "ERROR: invalid test (" << sodtest_num << ")." << endl;
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
  int rank, size; 
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  // set simulation parameters
  set_default_param(rank,size);
  parse_command_line_options(rank,size,argc,argv);
  set_param(rank,size); 

  // screen output
  if(rank==0){
    cout << "Sod test #" << sodtest_num << " in 1D:" << endl
         << " - number of particles: " << nparticles << endl
         << " - particles per core:  " << nparticlesproc << endl
         << " - output file: " << output_filename << endl;
  }

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

  // Header data 
  // the number of particles = nparticles 
  // The value for constant timestep 
  double timestep = 0.001;
  int dimension = 1;
  
  
  for(int64_t part=0; part<nparticlesproc; ++part){
    x[part] = lposition;

    if(x[part] > middle){
      P[part] = pressure_2;
      rho[part] = rho_2; 
      u[part] = u_2;
      //m[part] = m_2;
    }else{
      P[part] = pressure_1;
      rho[part] = rho_1;
      u[part] = u_1;
      //m[part] = m_1;
    }

    m[part] = rho[part]*middle/(nparticles/2.);

    //m[part] = 0.;
    // Y and Z not used 
    // VX, VY, VZ and AX, AY, AZ stay to 0
    h[part] = smoothing_length;
    
    // P stay to 0
    id[part] = posid++; 
    // Move to the next particle 
    lposition += ldistance;
    //std::cout<<x[part]<<": "<<h[part]<<std::endl;
  }

  // Destroy the file if exists 
  remove(output_filename);

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet; 
  testDataSet.createDataset(output_filename,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",nparticles);
  //testDataSet.writeDatasetAttribute("timestep","double",timestep);
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

  //_d1.createVariable("vx",Flecsi_Sim_IO::point,"double",nparticlesproc,vx);
  //_d2.createVariable("vy",Flecsi_Sim_IO::point,"double",nparticlesproc,vy);
  //_d3.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticlesproc,vz);

  //testDataSet.vars.push_back(_d1);
  //testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  //testDataSet.writeVariables();

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
  
  delete[] x,y,z,vx,vy,vz,ax,ay,az,h,rho,u,P,m,id,dt;
 
  MPI_Finalize();
  return 0;
}
