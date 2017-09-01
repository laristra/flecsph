/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 #include <iostream>
#include <algorithm>
#include <cmath>

#include "hdf5ParticleIO.h"


enum type {WALL=1, SIMPLE=0};

const double distance = 0.04;  // Distance between the particles 
const double m_ = 0.02;
const double rho_ = 998.29;
const double u_ = 1.0; 
const double smoothing_length = 0.0457;
const char* fileprefix = "hdf5_fluid_2D";
const double timestep = 0.01;
int32_t dimension = 2;

int main(int argc, char * argv[]){
  
  int rank, size; 
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  // Use only one process for generation in this version
  if(size > 1){
    if(rank==0){
      std::cerr<<
        "Use only one process for generation in this version"<<std::endl;
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  int nx = 10;//atoi(argv[1]);
  int ny = 10;//atoi(argv[2]);


  if(argc!=3){
    printf("./fluid_generator nx ny\n");
    fprintf(stderr,"Generation with default values= 10*10=100 particles");
  }else{
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
  }

  int64_t nparticles = nx*ny;
  int64_t nparticlesproc = nparticles/size;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  if(rank==0){
    printf("Generating %ld particles\n",nparticles);
    printf("%ld particles per proc (last %ld)\n",nparticlesproc,
        nparticles-nparticlesproc*(size-1));
  }
 
  // Generate data
  double nlines = nx; 
  double nlinesproc = nlines/size;
  double ncols = ny;
  double linestart = rank * nlinesproc * distance;

  double minX = -100.*distance;
  double minY = -1.*distance;
  double maxX = -2.*distance + nx*distance + 100.*distance; 
  double maxY = -2.*distance + ny*distance + 100.*distance;
  

  // Generate wall particles position 
  std::vector<double> vec_x; 
  std::vector<double> vec_y;
  //                    _  
  // Generate 4 lines | _ | for the tank 
  // Use two particles for wall like: 
  // o   o   o   o   o   
  //   o   o   o   o   o
  // Spaces by the smoothing length (not 2) 
    
  double div = 2.;
  // First left line 
  double val_x = minX - 10.*distance; 
  for(double i = minY-5.*distance; i < maxY+5.*distance ; i+=smoothing_length/div){
    vec_x.push_back(val_x);
    vec_y.push_back(i);  
    // Add the second particle 
    vec_x.push_back(val_x-sqrt(3./4.)*smoothing_length/div); 
    vec_y.push_back(i-1./2.*smoothing_length/div);
  }

  // Right line 
  val_x = maxX + 10.*distance; 
  for(double i = maxY+5.*distance ; i > minY-5.*distance ; i-=smoothing_length/div){
    vec_x.push_back(val_x);
    vec_y.push_back(i);  
    // Add the second particle 
    vec_x.push_back(val_x+sqrt(3./4.)*smoothing_length/div); 
    vec_y.push_back(i+1./2.*smoothing_length/div);
  }

  // Bottom line
  double val_y = minY - 5.*distance; 
  for(double i = maxX+10.*distance ; i>minX-10.*distance;i-=smoothing_length/div){
    vec_x.push_back(i);
    vec_y.push_back(val_y);  
    // Add the second particle 
    vec_x.push_back(i-1./2.*smoothing_length/div); 
    vec_y.push_back(val_y-sqrt(3./4.)*smoothing_length/div);
  }

  // Top line
  val_y = maxY + 5.*distance; 
  for(double i = minX-10.*distance; i<maxX+10.*distance;i+=smoothing_length/div){
    vec_x.push_back(i);
    vec_y.push_back(val_y);  
    // Add the second particle 
    vec_x.push_back(i+1./2.*smoothing_length/div); 
    vec_y.push_back(val_y+sqrt(3./4.)*smoothing_length/div);
  }

  int64_t nwall = vec_x.size(); 
  int64_t nwall_proc = nwall/size; 
  if(size==rank-1){
    nwall_proc = nwall - nwall_proc*(size-1);
  }

  if(rank==0){
    std::cout<<"Wall particles="<<nwall<<" per proc="<<nwall_proc <<std::endl; 
  }

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

  std::cout<<"Generating: "<<nlines<<"*"<<ncols<<std::endl;

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

    u[part] = u_; 
    rho[part] = rho_;
    m[part] = m_; 
    // Y and Z not used 
    // VX, VY, VZ and AX, AY, AZ stay to 0
    h[part] = smoothing_length;
    // P stay to 0
    id[part] = posid++; 

    type[part] = SIMPLE; 

    curcol += distance;
    col++;
    
    // Move to the next particle
    if(col == ncols){
      curline += distance; 
      curcol = 0.0;
      col = 0;
    }
  }

  // Then add the wall particles per process (this will shake the data) 
  int64_t pos = 0;
  for(int64_t i = (nwall/size)*rank ; i < (nwall/size)*rank+nwall_proc ; ++i){
    x[pos+nparticlesproc] = vec_x[i]; 
    y[pos+nparticlesproc] = vec_y[i]; 
    z[pos+nparticlesproc] = 0.;
    m[pos+nparticlesproc] = m_*100.;
    h[pos+nparticlesproc] = smoothing_length; 
    u[pos+nparticlesproc] = u_; 
    type[pos+nparticlesproc] = WALL; 
    rho[pos+nparticlesproc] = m[pos+nparticlesproc] / 
      (M_PI*4.*smoothing_length*smoothing_length);
    pos++;
  }

  char filename[128];
  sprintf(filename,"%s.h5part",fileprefix);
  std::cout<<"Writing to: "<<filename<<std::endl;
  remove(filename);

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(filename,MPI_COMM_WORLD);
  
  std::cout<<"Writing total of "<<nparticles+nwall<<" particles"<<std::endl
    <<std::flush;

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",nparticles+nwall);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",dimension);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",totalpart,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",totalpart,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",totalpart,z);

  std::cout<<"Pos writed"<<std::endl;

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",totalpart,vx);
  _d2.createVariable("vy",Flecsi_Sim_IO::point,"double",totalpart,vy);
  _d3.createVariable("vz",Flecsi_Sim_IO::point,"double",totalpart,vz);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("ax",Flecsi_Sim_IO::point,"double",totalpart,ax);
  _d2.createVariable("ay",Flecsi_Sim_IO::point,"double",totalpart,ay);
  _d3.createVariable("az",Flecsi_Sim_IO::point,"double",totalpart,az);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",totalpart,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",totalpart,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",totalpart,u);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",totalpart,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",totalpart,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",totalpart,id);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();
  
  _d1.createVariable("type",Flecsi_Sim_IO::point,"int32_t",totalpart,type);
    
  testDataSet.vars.push_back(_d1);
  
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
