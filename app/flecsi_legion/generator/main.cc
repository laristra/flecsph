#include <iostream>
#include <algorithm>

#include "hdf5ParticleIO.h"

const double distance = 0.003;  // Distance between the particles 
const double m_1 = 2.65e-3;
const double m_2 = 3.3125e-4;
const double rho_1 = 1;
const double rho_2 = 0.125;
const double u_1 = 0.0; //2.5;
const double u_2 = 0.0; //2;
const double smoothing_length = 1.0e-2;
const char* filename = "hdf5_sodtube.h5part";

int main(int argc, char * argv[]){

  if(argc!=2){
    printf("./sodtube_generator [nParticles]\n");
    exit(-1);
  }

  int rank, size; 
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int64_t nparticles = atoll(argv[1]);
  int64_t nparticlesproc = nparticles/size;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  if(rank==0){
    printf("Generating %ld particles\n",nparticles);
    printf("%ld particles per proc (last %ld)\n",nparticlesproc,
        nparticles-nparticlesproc*(size-1));
  }

  
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
  double middle = nparticles*distance/2.;
  // Find my first particle position 
  double position = distance*nparticlesproc*rank;
  // Id of my first particle 
  int64_t posid = nparticlesproc*rank;

  // Header data 
  // the number of particles = nparticles 
  // The value for constant timestep 
  double timestep = 0.001;
  int dimension = 1;
  
  
  for(int64_t part=0; part<nparticlesproc; ++part){
    x[part] = position;

    if(x[part] > middle){
      u[part] = u_2; 
      rho[part] = rho_2;
      m[part] = m_2; 
    }else{
      m[part] = m_1;
      u[part] = u_1;
      rho[part] = rho_1;
    }

    // Y and Z not used 
    // VX, VY, VZ and AX, AY, AZ stay to 0
    h[part] = smoothing_length;
    // P stay to 0
    id[part] = posid++; 
    // Move to the next particle 
    position += distance;
    //std::cout<<x[part]<<": "<<h[part]<<std::endl;
  }
 
 

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet; 
  testDataSet.createDataset(filename,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",nparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",1);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  char * simName = "sodtube_1D";
  testDataSet.writeDatasetAttributeArray("name","string",simName);
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
  _d2.createVariable("vy",Flecsi_Sim_IO::point,"double",nparticlesproc,vy);
  _d3.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticlesproc,vz);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("ax",Flecsi_Sim_IO::point,"double",nparticlesproc,ax);
  _d2.createVariable("ay",Flecsi_Sim_IO::point,"double",nparticlesproc,ay);
  _d3.createVariable("az",Flecsi_Sim_IO::point,"double",nparticlesproc,az);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();


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
