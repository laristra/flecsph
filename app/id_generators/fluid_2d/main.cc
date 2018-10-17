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

enum type {WALL=1, SIMPLE=0};

const double rest_density = 998.29;         // initial density
const double m_ = 0.02;                     // mass of an individual particle
const double rho_ = rest_density;
const double smoothing_length = 0.0457;
const double ldistance = smoothing_length;  // Distance between the particles 
const char* fileprefix = "hdf5_fluid_2D";
int32_t dimension = 2;

void generate_wall(
    std::vector<double>& x,
    std::vector<double>& y,
    std::vector<double>& vx, 
    std::vector<double>& vy,
    double x0, double y0, 
    double x1, double y1,
    double ldistance){
  double dist_points = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  clog(info)<<"Distance: "<<dist_points<<std::endl;
  int nelements = dist_points/ldistance;
  double xcur = x0;
  double ycur = y0;
  double xincr = (x1-x0)/(double)nelements;
  double yincr = (y1-y0)/(double)nelements;
  for(int i=0;i<nelements;++i){
    x.push_back(xcur);
    y.push_back(ycur);
    vx.push_back(0/*xincr*/);
    vy.push_back(0/*yincr*/);

    x.push_back(xcur+yincr);
    y.push_back(ycur-xincr);
    vx.push_back(0/*xincr*/);
    vy.push_back(0/*yincr*/);

    xcur += xincr; 
    ycur += yincr;
  }  
}

int main(int argc, char * argv[]){
  
  int rank, size;

  int provided; 
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  if(provided<MPI_THREAD_MULTIPLE){
	 clog(error)<<"Error MPI version: MPI_THREAD_MULTIPLE provided: "
    <<provided<<std::endl;
    MPI_Finalize(); 
	 exit(EXIT_FAILURE); 	
  }
 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // Use only one process for generation in this version
  if(size > 1){
    clog(error)<<
        "Use only one process for generation in this version"<<std::endl;
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  int nx = 10;//atoi(argv[1]);
  int ny = 10;//atoi(argv[2]);

  if(argc!=3){
    clog(warn)<<"./fluid_generator nx ny"<<std::endl;
    clog(warn)<<"Generation with default values= 10*10=100 particles"<<std::endl;
  }else{
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  int64_t nparticles = nx*ny;
  int64_t nparticlesproc = nparticles/size;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  clog(info)<<"Generating "<<nparticles<<" particles "<<std::endl;
  clog(info)<<nparticlesproc<<" particles per proc (last "<<
    nparticles-nparticlesproc*(size-1)<<")"<<std::endl;
 
  // Generate data
  double nlines = nx; 
  double nlinesproc = nlines/size;
  double ncols = ny;
  double linestart = rank * nlinesproc * ldistance;

  double minX = -2.*ldistance;
  double minY = -2.*ldistance;
  double maxX = nx*ldistance + 100.*ldistance; 
  double maxY = ny*ldistance + 100.*ldistance;
  

  // Generate wall particles position 
  std::vector<double> vec_x, vec_vx; 
  std::vector<double> vec_y, vec_vy;
  //                    _  
  // Generate 4 lines | _ | for the tank 
  // Use two particles for wall like: 
  // o   o   o   o   o   
  //   o   o   o   o   o
  // Spaces by the smoothing length (not 2) 
  
  //  |
  //  |
  //  |
  generate_wall(vec_x,vec_y,vec_vx,vec_vy,
      minX,minY-5.*ldistance,
      minX,maxY+5.*ldistance,
      smoothing_length); 
  //  |      |
  //  |      |
  //  |      |  
  generate_wall(vec_x,vec_y,vec_vx,vec_vy,
      maxX+10.*ldistance,minY-5.*ldistance,
      maxX+10.*ldistance,maxY+5.*ldistance,
      smoothing_length);
  //  |      |
  //  |      |
  //  |______|
  generate_wall(vec_x,vec_y,vec_vx,vec_vy,
      minX-10.*ldistance,minY,
      maxX+10.*ldistance,minY,
      smoothing_length);
  //  ________
  //  |      |
  //  |      |
  //  |______|  
  generate_wall(vec_x,vec_y,vec_vx,vec_vy,
      minX-10.*ldistance,maxY+5.*ldistance,
      maxX+10.*ldistance,maxY+5.*ldistance,
      smoothing_length);

  int64_t nwall = vec_x.size(); 
  int64_t nwall_proc = nwall/size; 
  if(size==rank-1){
    nwall_proc = nwall - nwall_proc*(size-1);
  }

  clog(info)<<"Wall particles="<<nwall<<" per proc="<<nwall_proc <<std::endl; 

  int64_t totalpart = nparticlesproc+nwall_proc; 

  // Position
  double* x = new double[totalpart]();
  double* y = new double[totalpart]();
  double* z = new double[totalpart]();
  // Velocity
  double* vx = new double[totalpart]();
  double* vy = new double[totalpart]();
  double* vz = new double[totalpart]();
  // Acceleration
  double* ax = new double[totalpart]();
  double* ay = new double[totalpart]();
  double* az = new double[totalpart]();
  // Smoothing length 
  double* h = new double[totalpart]();
  // Density 
  double* rho = new double[totalpart]();
  // Internal Energy 
  double* u = new double[totalpart]();
  // Pressure
  double* P = new double[totalpart]();
  // Mass
  double* m = new double[totalpart]();
  // Id
  int64_t* id = new int64_t[totalpart](); 
  // Type
  int32_t* type = new int32_t[totalpart]();

  // Timestep 
  double* dt = new double[totalpart]();

  clog(info)<<"Generating: "<<nlines<<"*"<<ncols<<std::endl;

  // Id of my first particle 
  int64_t posid = nparticlesproc*rank;

  // Header data 
  // the number of particles = nparticles 
  // The value for constant timestep 
   
  double curline = linestart;
  double curcol = 0.0;
  int col = 0;

  for(int64_t part=0; part<nparticlesproc; ++part){
    x[part] = curline;
    y[part] = curcol; 
    z[part] = 0.;

    rho[part] = rho_;
    m[part] = m_; 
    // Y and Z not used 
    // VX, VY, VZ and AX, AY, AZ stay to 0
    h[part] = smoothing_length;
    // P stay to 0
    id[part] = posid++; 

    type[part] = SIMPLE; 

    curcol += ldistance;
    col++;
    
    // Move to the next particle
    if(col == ncols){
      curline += ldistance; 
      curcol = 0.0;
      col = 0;
    }
  }

  // Then add the wall particles per process (this will shake the data) 
  int64_t pos = 0;
  for(int64_t i = (nwall/size)*rank ; i < (nwall/size)*rank+nwall_proc ; ++i){
    x[pos+nparticlesproc] = vec_x[i]; 
    y[pos+nparticlesproc] = vec_y[i]; 
    vx[pos+nparticlesproc] = vec_vx[i];
    vy[pos+nparticlesproc] = vec_vy[i];
    z[pos+nparticlesproc] = 0.;
    m[pos+nparticlesproc] = m_;
    h[pos+nparticlesproc] = smoothing_length; 
    type[pos+nparticlesproc] = WALL; 
    rho[pos+nparticlesproc] = rest_density;
    pos++;
  }

  char filename[128];
  sprintf(filename,"%s.h5part",fileprefix);
  clog(info)<<"Writing to: "<<filename<<std::endl;
  remove(filename);
  
  clog(info)<<"Writing total of "<<nparticles+nwall<<" particles"<<std::endl
    <<std::flush;

  h5_file_t * dataFile = H5OpenFile(filename,H5_O_WRONLY, MPI_COMM_WORLD);
    
  int use_fixed_timestep = 1; 
  // add the global attributes
  H5WriteFileAttribInt64(dataFile,"nparticles",&totalpart,1);
  double timestep = 0.001; 
  H5WriteFileAttribFloat64(dataFile,"timestep",&timestep,1);
  H5WriteFileAttribInt32(dataFile,"dimension",&dimension,1);
  H5WriteFileAttribInt32(dataFile,"use_fixed_timestep",&use_fixed_timestep,1);

  H5SetStep(dataFile,0);
  H5PartSetNumParticles(dataFile,totalpart);
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
