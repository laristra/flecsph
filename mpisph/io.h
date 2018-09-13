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
#include "params.h"

#include <H5hut.h>
#include <hdf5ParticleIO.h>

//#include "tree.h"
//#include "physics.h"

namespace io{

double 
input_parameter_double(
    const char * filename, 
    const char * attribute_name)
{

  auto dataFile = H5OpenFile(filename,H5_O_RDONLY
      | H5_VFD_MPIIO_IND,// Flag to be not use mpiposix
      MPI_COMM_WORLD);
  
  double value;
  if(H5_SUCCESS != H5ReadFileAttribFloat64(dataFile,attribute_name,&value)){
    value = double{};
  }
  H5CloseFile(dataFile);
  return value;
}


int 
input_parameter_int(
    const char * filename, 
    const char * attribute_name)
{

  auto dataFile = H5OpenFile(filename,H5_O_RDONLY
      | H5_VFD_MPIIO_IND,// Flag to be not use mpiposix
      MPI_COMM_WORLD);
  
  int value;
  if(H5_SUCCESS != H5ReadFileAttribInt32(dataFile,attribute_name,&value)){
    value = int{};
  }
  H5CloseFile(dataFile);
  return value;
}

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

  clog_one(trace)<<"Input particles";
  
  //--------------- OPEN FILE AND READ NUMBER OF PARTICLES -------------------- 
  auto dataFile = H5OpenFile(filename,H5_O_RDONLY
      | H5_VFD_MPIIO_IND,// Flag to be not use mpiposix
      MPI_COMM_WORLD);

  int val = H5HasStep(dataFile,startiteration);
  
  if(val != 1){
    clog_one(error)<<std::endl<<"Step "<<startiteration<<" not found in file "<<
        filename<<std::endl<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }
  // If yes, get into this step 
  H5SetStep(dataFile,startiteration);
  
  // Get the number of particles 
  int64_t nparticles = 0UL; 
  nparticles = H5PartGetNumParticles(dataFile);
  int64_t nparticlesproc = nparticles/size;
  // Handle the number of particles for the last one
  if(size-1==rank){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  // Register for main 
  totalnbodies = nparticles;
  nbodies = nparticlesproc;

  //--------------- READ GLOBAL ATTRIBUTES ------------------------------------
  // read the number of dimension 
  int32_t dimension;
  if(H5_SUCCESS == H5ReadFileAttribInt32(dataFile,"dimension",&dimension)){
    assert(gdimension == dimension);
  }else{
    clog_one(error)<<"No dimension value: setting to default 3"<<std::endl;
    dimension = 3;
  }

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
  int * dataInt32 = new int[nparticlesproc]; 

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

  // Reset buffer to 0, if next value not present 
  std::fill(dataX,dataX+nparticlesproc,0.);
  std::fill(dataY,dataY+nparticlesproc,0.);
  std::fill(dataZ,dataZ+nparticlesproc,0.); 

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
 
  // Reset buffer to 0, if next value not present 
  std::fill(dataX,dataX+nparticlesproc,0.);
  std::fill(dataY,dataY+nparticlesproc,0.);
  std::fill(dataZ,dataZ+nparticlesproc,0.); 

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
 

  // Reset buffer to 0, if next value not present 
  std::fill(dataX,dataX+nparticlesproc,0.);
  std::fill(dataY,dataY+nparticlesproc,0.);
  std::fill(dataZ,dataZ+nparticlesproc,0.); 

  H5PartReadDataFloat64(dataFile,"m",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setMass(dataX[i]);
  }
  H5PartReadDataFloat64(dataFile,"rho",dataY);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setDensity(dataY[i]);
  }
  H5PartReadDataFloat64(dataFile,"h",dataZ);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setSmoothinglength(dataZ[i]);
  }

  // Reset buffer to 0, if next value not present 
  std::fill(dataX,dataX+nparticlesproc,0.);
  std::fill(dataY,dataY+nparticlesproc,0.);
  std::fill(dataZ,dataZ+nparticlesproc,0.); 

  // Pressure  
  H5PartReadDataFloat64(dataFile,"P",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setPressure(dataX[i]);
  }

  // Internal Energy  
  #ifdef INTERNAL_ENERGY
  clog_one(trace)<<"Reading internal energy"<<std::endl;
  std::fill(dataX,dataX+nparticlesproc,0.);
  H5PartReadDataFloat64(dataFile,"u",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setInternalenergy(dataX[i]);
  }
  #endif

  //bool b_index = false; 
  // \TODO check if user ID is uniq
  bool b_index = H5_SUCCESS == 
  H5PartReadDataInt64(dataFile,"id",dataInt);
  if(b_index){
    for(int64_t i=0; i<nparticlesproc; ++i){
      bodies[i].second.setId(dataInt[i]); 
    }
  }else{
    clog_one(trace)<<"Setting ID for particles"<<std::endl;
    // Otherwise generate the id 
    int64_t start = (totalnbodies/size)*rank+1;
    for(int64_t i=0; i<nparticlesproc; ++i){
      bodies[i].second.setId(start+i); 
    }
    if(size-1 == rank){
      assert(totalnbodies==bodies.back().second.getId()); 
    }
  }
  
  // Reset buffer to 0, if next value not present 
  std::fill(dataX,dataX+nparticlesproc,0.);

  // delta T   
  H5PartReadDataFloat64(dataFile,"dt",dataX);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setDt(dataX[i]);
  }
  
  // Reset buffer to 0, if next value not present 
  std::fill(dataInt32,dataInt32+nparticlesproc,0.);

  // delta T   
  H5PartReadDataInt32(dataFile,"type",dataInt32);
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setType(dataInt32[i]);
  }


  delete[] dataX;
  delete[] dataY;
  delete[] dataZ;
  delete[] dataInt;
  delete[] dataInt32;

  MPI_Barrier(MPI_COMM_WORLD);
  H5CloseFile(dataFile);

  clog_one(trace)<<".done"<<std::endl;

}// inputDataHDF5

// Output data in HDF5 format 
// Generate the associate XDMF file 
void outputDataHDF5(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    const char* fileprefix,
    int step,
    double totaltime,
    bool do_diff_files = param::out_h5data_separate_iterations)
{
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  MPI_Barrier(MPI_COMM_WORLD);

  clog_one(trace)<<"Output particles"<<std::flush;

  int64_t nparticlesproc = bodies.size();

  char filename[128];
  if(do_diff_files){
    sprintf(filename,"%s_%05d.h5part",fileprefix,step*param::out_h5data_every);
  }
  else {
    sprintf(filename,"%s.h5part",fileprefix);
  }
  Flecsi_Sim_IO::HDF5ParticleIO simio;
  simio.createDataset(filename,MPI_COMM_WORLD);
  
  //-------------------GLOBAL HEADER-------------------------------------------
  // Only for the first output
  if(do_diff_files) {
     simio.writeDatasetAttribute("ndim","int32_t",gdimension);
     // ... TODO
  }
  else{
    if(step == 0){
      // output dimension 
      simio.writeDatasetAttribute("ndim","int32_t",gdimension);
    }
  }

  //------------------STEP HEADER----------------------------------------------
  // Put the step header
  simio.setTimeStep(step);
  
  Flecsi_Sim_IO::Attribute timeValue("time",Flecsi_Sim_IO::timestep,"double",
      totaltime);
  simio.addTimeStepAttribute(timeValue);


  simio.writeTimestepAttributes();

  //------------------STEP DATA------------------------------------------------

  // 4 buffers, 3 double and 1 int64
  double* b1 = new double[nparticlesproc];
  double* b2 = new double[nparticlesproc];
  double* b3 = new double[nparticlesproc];
  int64_t* bi = new int64_t[nparticlesproc];
  int32_t* bint = new int32_t[nparticlesproc];  


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

  // Smoothing length, Density, Internal Energy 
  pos = 0L;
  // Extract data from bodies 
  for(auto bi: bodies){
    b1[pos] = bi.second.getSmoothinglength();
    b2[pos] = bi.second.getDensity();
    #ifdef INTERNAL_ENERGY
    b3[pos] = bi.second.getInternalenergy();
    #endif
    pos++;
  }

  simio.addVariable( Flecsi_Sim_IO::Variable("h",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b1));
  simio.addVariable( Flecsi_Sim_IO::Variable("rho",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b2));
  #ifdef INTERNAL_ENERGY
  simio.addVariable( Flecsi_Sim_IO::Variable("u",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b3));
  #endif

  simio.writeVariables();

 // Pressure, Mass, Id, timestep
  pos = 0L;
  // Extract data from bodies 
  for(auto bid: bodies){
    b1[pos] = bid.second.getPressure();
    b2[pos] = bid.second.getMass();
    b3[pos] = bid.second.getDt();
    bi[pos] = bid.second.getId();
    bint[pos++] = bid.second.getType(); 
  }
  
  simio.addVariable( Flecsi_Sim_IO::Variable("P",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b1));
  simio.addVariable( Flecsi_Sim_IO::Variable("m",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b2));
  simio.addVariable( Flecsi_Sim_IO::Variable("dt",Flecsi_Sim_IO::point, 
        "double", nparticlesproc,b3));
  simio.addVariable( Flecsi_Sim_IO::Variable("id",Flecsi_Sim_IO::point, 
        "int64_t", nparticlesproc,bi));
  simio.addVariable( Flecsi_Sim_IO::Variable("type",Flecsi_Sim_IO::point, 
        "int32_t", nparticlesproc,bint));


  simio.writeVariables();
  simio.closeFile();

  delete[] b1;
  delete[] b2;
  delete[] b3;
  delete[] bi;
  delete[] bint;

  clog_one(trace)<<".done"<<std::endl;
  
}// outputDataHDF5

} // namespace io

#endif


