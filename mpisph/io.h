/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file io.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Simple implementation for Input Output for serial and distributed
 */

#ifndef _mpisph_io_h_
#define _mpisph_io_h_

#include <cstdlib>
#include <iostream>
#include <vector>

#include <H5hut.h>
#include <hdf5ParticleIO.h>

//#include "tree.h"
#include "physics/physics.h"

namespace io{

// Read data from txt file with ranges
// This is not the most beautiful way, but enough for testing
void inputDataTxtRange(
    std::vector<std::pair<entity_key_t,body>>& bodies, 
    int64_t& nbodies,
    int64_t& totalnbodies,
    const char * filename
  )
{
  int rank, size; 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  nbodies = 0; 
  totalnbodies = 0;
  int mpi_err;
  MPI_File file; 
  MPI_Offset filesize;
  MPI_Offset nelemperproc;
  MPI_Offset startposition; 
  MPI_Offset stopposition; 
  double x, y, z, velX, velY, velZ, velHX, velHY, velHZ;
  double accX, accY, accZ, smoothinglength, pressure;
  double entropy, density, mass, tmp, angularMoment;
  int tmpD;

  //int overlap = 100;

  mpi_err = MPI_File_open(MPI_COMM_WORLD,filename,
      MPI_MODE_RDONLY,MPI_INFO_NULL, &file);
  if(mpi_err){
    fprintf(stderr,"Could not open file %s\n",filename);
    MPI_Finalize(); 
    exit(1); 
  }

  MPI_File_get_size(file,&filesize);
  nelemperproc = filesize/size;
  startposition = rank * nelemperproc; 
  stopposition = startposition + nelemperproc-1;
  if(rank == size-1)
    stopposition = filesize;
  //else 
  //  stopposition += overlap; 
  // Read in a buffer of a double and or int 
  char buffer[2048]; 
  //int iintbuffer; 
 
  //printf("%d/%d from %d to %d",rank,size,startposition,stopposition);

  // Read a line in char buffer
  MPI_File_read_at(file,startposition,
      &buffer,2048,MPI_CHAR,MPI_STATUS_IGNORE); 
  // Print the read value 
  //std::cout << "Just read: "<< buffer<< std::endl;
  
  // For all but the first process, go to the end of this line 
  if(rank != 0){
    int pos = 0;
    while(buffer[pos]!='\n'){
      startposition++;
      pos++;
    }
    startposition++;
    // Read another line from new starting point
    MPI_File_read_at(file,startposition,
      &buffer,2048,MPI_CHAR,MPI_STATUS_IGNORE); 
    // Print the read value 
    //std::cout << "Second read: "<< buffer<< std::endl;
  }

  // Interpret the buffer 
  sscanf(buffer,
      "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
      " %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf \n",
      &x,&y,&z,                       // positions
      &velX,&velY,&velZ,              // velocity
      &velHX,&velHY,&velHZ,           // velHalf 
      &accX,&accY,&accZ,              // acceleration
      &density,&pressure,&entropy,    // density, pressure, entropy 
      &mass,&smoothinglength,&tmp,    // mass, smoothing length, dt 
      &tmpD,&tmp,&tmp,                // step, totalTime, tmp
      &tmp,&tmp,&angularMoment        // tmp, tmp, angularMoment
  );
  point_t position = {x,y,z}; 
  point_t velocity = {velX,velY,velZ};
  point_t velocityhalf = {velHX,velHY,velHZ};
  point_t acceleration = {accX, accY, accZ};
  
  auto bi = body(position,velocity,velocityhalf,
      acceleration,density,pressure,entropy,mass,smoothinglength);
  bodies.push_back(std::make_pair(entity_key_t::null(),bi));
  ++nbodies;

  // Move at the end of this line 
  int pos = 0;
  while(buffer[pos]!='\n'){
    startposition++;
    pos++;
  }
  startposition++;

  // Loop for all the lines of the file 
  while(startposition < stopposition){
    // Read new line  
    MPI_File_read_at(file,startposition,
      &buffer,2048,MPI_CHAR,MPI_STATUS_IGNORE); 
    // Interpret
    sscanf(buffer,
      "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
      " %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf \n",
      &x,&y,&z,                       // positions
      &velX,&velY,&velZ,              // velocity
      &velHX,&velHY,&velHZ,           // velHalf 
      &accX,&accY,&accZ,              // acceleration
      &density,&pressure,&entropy,    // density, pressure, entropy 
      &mass,&smoothinglength,&tmp,    // mass, smoothing length, dt 
      &tmpD,&tmp,&tmp,                // step, totalTime, tmp
      &tmp,&tmp,&angularMoment        // tmp, tmp, angularMoment
    );
    position = {x,y,z}; 
    velocity = {velX,velY,velZ};
    velocityhalf = {velHX,velHY,velHZ};
    acceleration = {accX, accY, accZ};

    auto bi = body(position,velocity,velocityhalf,
      acceleration,density,pressure,entropy,mass,smoothinglength);
    bodies.push_back(std::make_pair(entity_key_t::null(),bi));
    
    ++nbodies;
    //std::cout << *bi << std::endl;
    // Move at the end of this line 
    int pos = 0;
    while(buffer[pos]!='\n'){
      startposition++;
      pos++;
    }
    startposition++;
  }
  // Reduction over the bodies 
  MPI_Allreduce(&nbodies,&totalnbodies,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  std::cout << rank<<"/"<<size<<" nbodies: "<<nbodies<<"/"<<totalnbodies<<std::endl; 

} // inputDataTxtRange

// Input data fro HDF5 File
void inputDataHDF5(
  std::vector<std::pair<entity_key_t,body>>& bodies,
  const char * filename,
  int64_t& totalnbodies,
  int64_t& nbodies,
  int startiteration)
{
  
  int rank, size; 
  MPI_Comm_size(MPI_COMM_WORLD,&size); 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    std::cout<<"Input particles";
  }
  MPI_Barrier(MPI_COMM_WORLD); 

  //--------------- OPEN FILE AND READ NUMBER OF PARTICLES -------------------- 
  auto dataFile = H5OpenFile(filename,H5_O_RDONLY|
      H5_VFD_INDEPENDENT, // Flag to be not use mpiposix
      MPI_COMM_WORLD);
  // Check if timestep exists
  int val = H5HasStep(dataFile,startiteration);
  if(val != 1){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      std::cout<<std::endl<<"Step "<<startiteration<<" not found in file "<<
        filename<<std::endl<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }
  // If yes, get into this step 
  H5SetStep(dataFile,startiteration);
   
  // Get the number of particles 
  int64_t nparticles = 0UL; 
  nparticles = H5PartGetNumParticles(dataFile);
  if(rank == 0){
    std::cout<<"Total of "<<nparticles<<" particles"<<std::endl;
  }
  int64_t nparticlesproc = nparticles/size;
  // Handle the number of particles for the last one
  if(size-1==rank){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  // Register for main 
  totalnbodies = nparticles;
  nbodies = nparticlesproc;

  std::cout<<rank<<": "<<nparticlesproc<<std::endl;

  //--------------- READ GLOBAL ATTRIBUTES ------------------------------------
  // read the number of dimension 
  int32_t dimension;
  H5ReadFileAttribInt32(dataFile,"dimension",&dimension);
  assert(gdimension == dimension);

  // timestep 
  H5ReadFileAttribFloat64(dataFile,"timestep",&physics::dt);

  //--------------- READ DATA FROM STEP ---------------------------------------
  // Set the number of particles read by each process
  H5PartSetNumParticles(dataFile,nparticlesproc);
  // Resize the bodies vector 
  bodies.clear();
  bodies.resize(nparticlesproc); 

  // Read the dataset and fill the particles data 
  double* dataX = new double[nparticlesproc];
  double* dataY = new double[nparticlesproc];
  double* dataZ = new double[nparticlesproc];
  int64_t* dataInt = new int64_t[nparticlesproc];

  std::cout<<rank<<": array allocated"<<std::endl;

  // Positions
  H5PartReadDataFloat64(dataFile,"x",dataX);
  H5PartReadDataFloat64(dataFile,"y",dataY);
  H5PartReadDataFloat64(dataFile,"z",dataZ);
  for(int64_t i=0; i<nparticlesproc; ++i){
    point_t position; 
    position[0] = dataX[i];
    if(gdimension>1){
      position[1] = dataY[i];
    }
    if(gdimension>2){
      position[2] = dataZ[i];
    }
    bodies[i].second.setPosition(position);
  }

  std::cout<<rank<<": position read"<<std::endl;

  // Velocity
  H5PartReadDataFloat64(dataFile,"vx",dataX);
  H5PartReadDataFloat64(dataFile,"vy",dataY);
  H5PartReadDataFloat64(dataFile,"vz",dataZ);
  for(int64_t i=0; i<nparticlesproc; ++i){
    point_t velocity;
    velocity[0] = dataX[i];
    if(gdimension>1){
      velocity[1] = dataY[i];
    }
    if(gdimension>2){
      velocity[2] = dataZ[i];
    }
    bodies[i].second.setVelocity(velocity);
  }

  // Acceleration
  H5PartReadDataFloat64(dataFile,"ax",dataX);
  H5PartReadDataFloat64(dataFile,"ay",dataY);
  H5PartReadDataFloat64(dataFile,"az",dataZ);
  for(int64_t i=0; i<nparticlesproc; ++i){
    point_t acceleration;
    acceleration[0] = dataX[i];
    if(gdimension>1){
      acceleration[1] = dataY[i];
    }
    if(gdimension>2){
      acceleration[2] = dataZ[i];
    }
    bodies[i].second.setAcceleration(acceleration);
  }

  // Smoothing Length 
  H5PartReadDataFloat64(dataFile,"h",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setSmoothinglength(dataX[i]);
  }

  // Density   
  H5PartReadDataFloat64(dataFile,"rho",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setDensity(dataX[i]);
  }

  // Internal Energy   
  H5PartReadDataFloat64(dataFile,"u",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setInternalenergy(dataX[i]);
  }

  // Pressure  
  H5PartReadDataFloat64(dataFile,"P",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setPressure(dataX[i]);
  }
 
 // Mass  
  H5PartReadDataFloat64(dataFile,"m",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setMass(dataX[i]);
  }

  // Id   
  H5PartReadDataInt64(dataFile,"id",dataInt);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setId(dataInt[i]);
  }
 
  // delta T   
  H5PartReadDataFloat64(dataFile,"dt",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setDt(dataX[i]);
  }

  std::cout<<rank<<": done last"<<std::endl;
  MPI_Barrier(MPI_COMM_WORLD);


  delete[] dataX;
  delete[] dataY;
  delete[] dataZ;
  delete[] dataInt;

  MPI_Barrier(MPI_COMM_WORLD);
  H5CloseFile(dataFile);

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    std::cout<<".done"<<std::endl;
  }

}// inputDataHDF5

// Output data in HDF5 format 
// Generate the associate XDMF file 
void outputDataHDF5(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    const char* fileprefix,
    int step,
    bool do_diff_files = false)
{

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0){
    std::cout<<"Output particles";
  }

  int64_t nparticlesproc = bodies.size();

  char filename[128];
  if(do_diff_files){
    sprintf(filename,"%s_%d.h5part",fileprefix,step);
  }else{
    sprintf(filename,"%s.h5part",fileprefix);
  }  

  // open the file 
  // auto dataFile = H5OpenFile(filename,H5_O_RDWR,MPI_COMM_WORLD);
  // Set the step 
  // H5SetStep(dataFile,step);
  // Set the number of particles 
  // H5PartSetNumParticles(dataFile,nparticlesproc);

  Flecsi_Sim_IO::HDF5ParticleIO simio;
  simio.createDataset(filename,MPI_COMM_WORLD);

  //-------------------GLOBAL HEADER-------------------------------------------
  // Only for the first output

  //------------------STEP HEADER----------------------------------------------
  // Put the step header
  simio.setTimeStep(step);
  simio.addTimeStepAttribute(
      Flecsi_Sim_IO::Attribute(
        "dt",
        Flecsi_Sim_IO::timestep,
        "double",
        physics::totaltime)
      );


  //------------------STEP DATA------------------------------------------------

  // 4 buffers, 3 double and 1 int64
  double* b1 = new double[nparticlesproc];
  double* b2 = new double[nparticlesproc];
  double* b3 = new double[nparticlesproc];
  int64_t* bi = new int64_t[nparticlesproc];

  // Position
  int64_t pos = 0L;
  // Extract data from bodies 
  for(auto bi: bodies){
    b1[pos] = bi.second.getPosition()[0];
    if(gdimension>1){
      b2[pos] = bi.second.getPosition()[1];
    }else{
      b2[pos] = 0.;
    }
    if(gdimension>2){
      b3[pos++] = bi.second.getPosition()[2];
    }else{
      b3[pos++] = 0.;
    }
  }

  // Add variable  
  simio.addVariable( Flecsi_Sim_IO::Variable("x",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b1));
  simio.addVariable( Flecsi_Sim_IO::Variable("y",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b2));
  simio.addVariable( Flecsi_Sim_IO::Variable("z",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b3));
  // Push to file 
  simio.writeVariables();

  // Velocity
  pos = 0L;
  // Extract data from bodies 
  for(auto bi: bodies){
    b1[pos] = bi.second.getVelocity()[0];
    if(gdimension>1){
      b2[pos] = bi.second.getVelocity()[1];
    }else{
      b2[pos] = 0;
    }
    if(gdimension>2){
      b3[pos++] = bi.second.getVelocity()[2];
    }else{
      b3[pos++] = 0.;
    }
  }

  simio.addVariable( Flecsi_Sim_IO::Variable("vx",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b1));
  simio.addVariable( Flecsi_Sim_IO::Variable("vy",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b2));
  simio.addVariable( Flecsi_Sim_IO::Variable("vz",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b3));

  //H5PartWriteDataFloat64(dataFile,"vx",dataX);
  //H5PartWriteDataFloat64(dataFile,"vy",dataY);
  //H5PartWriteDataFloat64(dataFile,"vz",dataZ);

  simio.writeVariables();

  // Acceleration 
  pos = 0L;
  // Extract data from bodies 
  for(auto bi: bodies){
    b1[pos] = bi.second.getAcceleration()[0];
    if(gdimension>1){
      b2[pos] = bi.second.getAcceleration()[1];
    }else{
      b2[pos] = 0.;
    }
    if(gdimension>2){
      b3[pos++] = bi.second.getAcceleration()[2]; 
    }else{
      b3[pos++] = 0.;
    }
  }

  simio.addVariable( Flecsi_Sim_IO::Variable("ax",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b1));
  simio.addVariable( Flecsi_Sim_IO::Variable("ay",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b2));
  simio.addVariable( Flecsi_Sim_IO::Variable("az",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b3));

  simio.writeVariables();

  // H5PartWriteDataFloat64(dataFile,"vx",dataX);
  // H5PartWriteDataFloat64(dataFile,"vy",dataY);
  // H5PartWriteDataFloat64(dataFile,"vz",dataZ);

  // Smoothing length, Density, Internal Energy 
  pos = 0L;
  // Extract data from bodies 
  for(auto bi: bodies){
    b1[pos] = bi.second.getSmoothinglength();
    b2[pos] = bi.second.getDensity();
    b3[pos++] = bi.second.getInternalenergy();
  }

  simio.addVariable( Flecsi_Sim_IO::Variable("h",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b1));
  simio.addVariable( Flecsi_Sim_IO::Variable("rho",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b2));
  simio.addVariable( Flecsi_Sim_IO::Variable("u",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b3));

  simio.writeVariables();

  // H5PartWriteDataFloat64(dataFile,"h",dataX);
  // H5PartWriteDataFloat64(dataFile,"rho",dataY);
  // H5PartWriteDataFloat64(dataFile,"u",dataZ);

 // Pressure, Mass, Id, timestep
  pos = 0L;
  // Extract data from bodies 
  for(auto bid: bodies){
    b1[pos] = bid.second.getPressure();
    b2[pos] = bid.second.getMass();
    b3[pos] = bid.second.getDt();
    bi[pos++] = bid.second.getId();
  }
  
  simio.addVariable( Flecsi_Sim_IO::Variable("P",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b1));
  simio.addVariable( Flecsi_Sim_IO::Variable("m",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b2));
  simio.addVariable( Flecsi_Sim_IO::Variable("dt",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b3));
  simio.addVariable( Flecsi_Sim_IO::Variable("id",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,bi));

  simio.writeVariables();
  simio.closeFile();
  // H5PartWriteDataFloat64(dataFile,"P",dataX);
  // H5PartWriteDataFloat64(dataFile,"m",dataY);
  // H5PartWriteDataFloat64(dataFile,"dt",dataZ);
  // H5PartWriteDataInt64(dataFile,"id",dataInt);

  //int64_t nparticles;
  //MPI_Allreduce(&nparticles,&nparticlesproc,1,MPI_INT64_T,MPI_COMM_WORLD);

  // Also output nparticles and timestep 
  //H5PartWriteDataInt64(dataFile,"nparticles",&nparticles);

  // H5CloseFile(dataFile);

  delete[] b1;
  delete[] b2;
  delete[] b3;
  delete[] bi;

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0){
    std::cout<<".done"<<std::endl;
  }

}// outputDataHDF5

} // namespace io

#endif


