#include <iostream>
#include <algorithm>

#include "hdf5SimIO.h"

const double distance = 0.0026566;  // Distance between the particles 
const double m_1 = 2.65e-3;
const double m_2 = 3.3125e-4;
const double rho_1 = 1;
const double rho_2 = 0.125;
const double u_1 = 2.5;
const double u_2 = 2;
const double smoothing_length = 1.0e-2;
const char* filename = "hdf5_sodtube.h5";

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
 
  
  Flecsi_Sim_IO::HDF5SimIO testDataSet; 
  testDataSet.setFilename(filename);
  testDataSet.createDataset(Flecsi_Sim_IO::unStructuredGrid,18,3);

  Flecsi_Sim_IO::Variable _x,_y,_z,_vx,_vy,_vz,_ax,_ay,_az,
    _h,_rho,_u,_P,_m,_id,_dt, _nparticles, _timestep;

  _x.createVariable("x",Flecsi_Sim_IO::point,"double",nparticlesproc,x);
  _y.createVariable("y",Flecsi_Sim_IO::point,"double",nparticlesproc,y);
  _z.createVariable("z",Flecsi_Sim_IO::point,"double",nparticlesproc,z);

  _vx.createVariable("vx",Flecsi_Sim_IO::point,"double",nparticlesproc,vx);
  _vy.createVariable("vy",Flecsi_Sim_IO::point,"double",nparticlesproc,vy);
  _vz.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticlesproc,vz);

  _ax.createVariable("ax",Flecsi_Sim_IO::point,"double",nparticlesproc,ax);
  _ay.createVariable("ay",Flecsi_Sim_IO::point,"double",nparticlesproc,ay);
  _az.createVariable("az",Flecsi_Sim_IO::point,"double",nparticlesproc,az);

  _h.createVariable("h",Flecsi_Sim_IO::point,"double",nparticlesproc,h);
  _rho.createVariable("rho",Flecsi_Sim_IO::point,"double",nparticlesproc,rho);
  _u.createVariable("u",Flecsi_Sim_IO::point,"double",nparticlesproc,u);
  _P.createVariable("P",Flecsi_Sim_IO::point,"double",nparticlesproc,P);
  _m.createVariable("m",Flecsi_Sim_IO::point,"double",nparticlesproc,m);
  _id.createVariable("id",Flecsi_Sim_IO::point,"int64_t",nparticlesproc,id);
  _dt.createVariable("dt",Flecsi_Sim_IO::point,"double",nparticlesproc,dt);
 
  _nparticles.createVariable("nparticles",
      Flecsi_Sim_IO::point,"int64_t",1,&nparticles);
  _timestep.createVariable("timestep",
      Flecsi_Sim_IO::point,"double",1,&timestep);

  testDataSet.vars.push_back(_x);
  testDataSet.vars.push_back(_y);
  testDataSet.vars.push_back(_z);

  testDataSet.vars.push_back(_vx);
  testDataSet.vars.push_back(_vy);
  testDataSet.vars.push_back(_vz);
  
  testDataSet.vars.push_back(_ax);
  testDataSet.vars.push_back(_ay);
  testDataSet.vars.push_back(_az);
  
  testDataSet.vars.push_back(_h);
  testDataSet.vars.push_back(_rho);
  testDataSet.vars.push_back(_u);
  testDataSet.vars.push_back(_P);
  testDataSet.vars.push_back(_m);
  testDataSet.vars.push_back(_id);
  testDataSet.vars.push_back(_dt);
  
  testDataSet.vars.push_back(_nparticles);
  testDataSet.vars.push_back(_timestep);

  testDataSet.writePointData(0,MPI_COMM_WORLD);

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
