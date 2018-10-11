/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <H5hut.h>


#include "user.h"
#include "fluid.h"
#include "kernels.h"


const double ldistance = 0.05;  // Distance between the particles 
const double m_ = 1.0e-5;
const double rho_ = 0.0;
const double u_ = 2.0; 
const double smoothing_length = 5e-2;
const char* fileprefix = "hdf5_fluid";
const double timestep = 1e-3;
int32_t dimension = 3;

int main(int argc, char * argv[]){
 
  int rank, size; 
  int provided; 
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE); 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  int nx = 10;//atoi(argv[1]);
  int ny = 10;//atoi(argv[2]);
  int nz = 10;//atoi(argv[3]);


  if(argc!=4){
    clog(warn)<<"./fluid_generator nx ny nz"<<std::endl;
    clog(warn)<<"Generation with default values= 10*10*10=1000 particles"
      <<std::endl;
  }else{
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nz = atoi(argv[3]);
  }

  int64_t nparticles = nx*ny*nz;
  int64_t nparticlesproc = nparticles/size;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  clog(info)<<"Generating "<<nparticles<<" particles"<<std::endl;
  clog(info)<<nparticlesproc<<" particles per proc (last "<<
    nparticles-nparticlesproc*(size-1)<<")"<<std::endl;

  if(nz == 0){
    dimension = 2;
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
  double nlines = nx; 
  double nlinesproc = nlines/size;
  double ncols = ny;
  double linestart = rank * nlinesproc * ldistance;
  double ndepth = nz;

  clog(info)<<"Generating: "<<nlines<<"*"<<ncols<<std::endl;

  // Id of my first particle 
  int64_t posid = nparticlesproc*rank;

  // Header data 
  // the number of particles = nparticles 
  // The value for constant timestep 
   
  double curline = linestart;
  double curcol = 0.0;
  double curdepth = 0.0;
  int col = 0;
  int depth = 0;

  for(int64_t part=0; part<nparticlesproc; ++part){
    x[part] = curline;
    y[part] = curcol; 
    z[part] = curdepth;

    u[part] = u_; 
    rho[part] = rho_;
    m[part] = m_; 
    // Y and Z not used 
    // VX, VY, VZ and AX, AY, AZ stay to 0
    h[part] = smoothing_length;
    // P stay to 0
    id[part] = posid++; 
    
    curcol += ldistance;
    //curdepth += distance; 
    col++;
    //depth++;
    // Move to the next particle
    if(col == ncols){
      curdepth += ldistance; 
      depth++;
      col = 0;
      curcol = 0.0;
      if( depth == ndepth  ){
        curline += ldistance; 
        curcol = 0.0;
        curdepth = 0.0;
        col = 0;
        depth = 0;
      }
    }
  }
  
  char filename[128];
  sprintf(filename,"%s.h5part",fileprefix);
  remove(filename);

  h5_file_t * dataFile = H5OpenFile(filename,H5_O_WRONLY, MPI_COMM_WORLD);
    
  int use_fixed_timestep = 1; 
  // add the global attributes
  H5WriteFileAttribInt64(dataFile,"nparticles",&nparticlesproc,1);
  H5WriteFileAttribFloat64(dataFile,"timestep",&timestep,1);
  H5WriteFileAttribInt32(dataFile,"dimension",&dimension,1);
  H5WriteFileAttribInt32(dataFile,"use_fixed_timestep",&use_fixed_timestep,1);

  H5SetStep(dataFile,0);
  H5PartSetNumParticles(dataFile,nparticlesproc);
  H5PartWriteDataFloat64(dataFile,"x",x);
  H5PartWriteDataFloat64(dataFile,"y",y);
  H5PartWriteDataFloat64(dataFile,"z",z);
  H5PartWriteDataFloat64(dataFile,"vx",vx);
  H5PartWriteDataFloat64(dataFile,"vy",vy);
  H5PartWriteDataFloat64(dataFile,"vz",vz);
  H5PartWriteDataFloat64(dataFile,"ax",ax);
  H5PartWriteDataFloat64(dataFile,"ay",ay);
  H5PartWriteDataFloat64(dataFile,"az",az);
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
