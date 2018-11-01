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
#include "default_physics.h"

#include <H5hut.h>

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
  const char * output_file_prefix, 
  int64_t& totalnbodies,
  int64_t& nbodies,
  int startIteration)
{

  char output_filename[128];
  sprintf(output_filename,"%s.h5part",output_file_prefix);

  // Default if new file, startStep = 0
  int startStep = 0;

  // Lower message level, handle them here on one process 
  H5SetVerbosityLevel  (0); 

  int rank, size; 
  MPI_Comm_size(MPI_COMM_WORLD,&size); 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  rank|| clog(trace)<<"Input particles" << std::endl;
 
  // The input file depend of the type of output, separate or one file. 
  h5_file_t * dataFile = nullptr; 
  if(startIteration == 0){
    dataFile = H5OpenFile(filename,H5_O_RDONLY
      | H5_VFD_MPIIO_IND,// Flag to be not use mpiposix
      MPI_COMM_WORLD); 
  }

  // ------------- CHECK IF THE SEARCHED ITERATION EXISTS ---------------------
  // Go through all the steps of the file and try to read the iteration 
  if(!param::out_h5data_separate_iterations && startIteration != 0){
   dataFile = H5OpenFile(filename,H5_O_RDONLY
      | H5_VFD_MPIIO_IND,// Flag to be not use mpiposix
      MPI_COMM_WORLD);
   int step = 1; 
    bool end = false; 
    bool found = false; 
    while(!end){
      int hasStep = H5HasStep(dataFile,step); 
      if(hasStep){
        // Check iteration number 
        H5SetStep(dataFile,step); 
        int64_t iteration; 
        if(H5_SUCCESS != H5ReadStepAttribInt64(dataFile,"iteration",
            &iteration)){
          rank || clog(error) << "Cannot read iteration in step "
            <<step<<std::endl;
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Finalize(); 
        }
        if(iteration == startIteration){
          found = true;
          end = true;  
        }
      }else{
        end = true; 
      }
      ++step; 
    }
    if(!found){
      rank || clog(error) << "Cannot find iteration "<<startIteration<<" in "
        <<filename<<std::endl;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize(); 
    }
    startStep = step-1;
    rank || clog(warn)<<"Found step "<<startStep<<" for iteration "<<
      startIteration<<std::endl; 
  }


  // Check if the file exists in case of multiple files
  if(param::out_h5data_separate_iterations && startIteration != 0){
    char step_filename[128];
    bool end = false; 
    bool found = false; 
    int step = 1;  
    while(!end){
      // Generate the filename associate with this step 
      sprintf(step_filename,"%s_%05d.h5part",output_file_prefix,step);
      rank || clog(trace) <<"Checking if file "<<step_filename<<" exists"
        <<std::endl<<std::flush;
      MPI_Barrier(MPI_COMM_WORLD); 
      // Check if files exists 
      if(access(step_filename,F_OK)==-1){
        rank || clog(error)<<"Cannot find file "<< step_filename<<
          " unable to find file with iteration "<<startIteration<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize(); 
      }
      // File exists, check the iteration 
      auto stepFile = H5OpenFile(step_filename,H5_O_RDONLY
          | H5_VFD_MPIIO_IND, MPI_COMM_WORLD);
      int64_t iteration; 
      H5SetStep(stepFile,step); 
      if(H5_SUCCESS != H5ReadStepAttribInt64(stepFile,"iteration",
            &iteration)){
          rank || clog(error) << "Cannot read iteration in file "
            <<step_filename<<std::endl;
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Finalize(); 
      }
      if(iteration == startIteration){
        found = true; 
        end = true; 
      }
      ++step; 
    }
    startStep = step-1;
    char diff_filename[128];
    sprintf(diff_filename,"%s_%05d.h5part",output_file_prefix,startStep);
    rank || clog(warn) << "Reading from file "<<diff_filename<<std::endl;
    // Change input file in this case 
    dataFile = H5OpenFile(diff_filename,H5_O_RDONLY
      | H5_VFD_MPIIO_IND,// Flag to be not use mpiposix
      MPI_COMM_WORLD);
  }

  // ------------- CHECK THAT THE STEP FOUND IS THE LAST IN OUTPUT ------------
  // If we are in single file mode 
  if(!param::out_h5data_separate_iterations){
    // If I start a new simulation, just delete the file
    // if not equal to the input file 
    if(startIteration == 0 && (strcmp(output_filename,filename) != 0)){
      remove(output_filename);
    }
    // If the file exists (either same as input or different)
    // Check if the lastStep is the startIteration 
    if(access(output_filename,F_OK)!=-1){
      auto outputFile = H5OpenFile(output_filename,H5_O_RDONLY
          | H5_VFD_MPIIO_IND, MPI_COMM_WORLD);
      // Check if startStep == lastStep 
      int lastStep = startStep;
      while(H5_SUCCESS == H5HasStep(outputFile,lastStep++)){}
      --lastStep; 
      rank || clog(trace)<<"startStep: "<<startStep<<" lastStep: "
        <<lastStep<<std::endl;
      if(startStep != lastStep){
        rank || clog(error) << "First step not last step in output"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize(); 
      }
    }
  }

  if(dataFile == nullptr){
    rank || clog(error) << "Cannot find data file"<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); 
  }

  int val = H5HasStep(dataFile,startStep);
  
  if(val != 1){
    rank|| clog(error)<<"Step "<<startStep<<" not found in file "<<
        filename<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }
 
  // Check if there is one step ahead of this one
  int step_ahead = H5HasStep(dataFile,startStep+1); 
  if(step_ahead == 1)
  {
    rank || clog(error)<<"ERROR Input file already have"<<
      " next step: "<<startStep+1<<std::endl<<std::flush;
    MPI_Barrier(MPI_COMM_WORLD); 
    MPI_Finalize(); 
  }

  // If yes, get into this step 
  H5SetStep(dataFile,startStep);
  
  // Get the number of particles 
  int64_t nparticles = 0UL; 
  nparticles = H5PartGetNumParticles(dataFile);
  int64_t nparticlesproc = nparticles/size;
  // Handle the number of particles for the last one
  if(size == rank + 1){
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
    rank|| clog(error)<<"No dimension value: setting to default 3"<<std::endl;
    dimension = 3;
  }

  //--------------- READ DATA FROM STEP ---------------------------------------
  // Set the number of particles read by each process
  H5PartSetNumParticles(dataFile,nparticlesproc);
  // Resize the bodies vector 
  bodies.clear();
  bodies.resize(nparticlesproc); 

  // Try to read timestep if exists 
  if(startIteration != 0){
    double timestep = 0.;
    double totaltime = 0.;
    if(H5_SUCCESS == H5ReadStepAttribFloat64(dataFile,"timestep",&timestep)){
      physics::dt = timestep; 
    }else{
      rank || clog(warn)<<"Unable to read timestep from file"<<std::endl;
    }
    if(H5_SUCCESS == H5ReadStepAttribFloat64(dataFile,"time",&totaltime)){
      physics::totaltime = totaltime; 
    }else{
      rank || clog(warn)<<"Unable to read totaltime from file"<<std::endl;
    }
  }

  // Read the dataset and fill the particles data 
  double* dataX = new double[nparticlesproc];
  double* dataY = new double[nparticlesproc];
  double* dataZ = new double[nparticlesproc];
  int64_t* dataInt = new int64_t[nparticlesproc];
  int * dataInt32 = new int[nparticlesproc];

  // Handle errors from H5HUT 
  h5_err_t errX, errY, errZ;  
  errX = errY = errZ = H5_SUCCESS; // prevent warnings

  // Positions
  errX = H5PartReadDataFloat64(dataFile,"x",dataX);
  if(gdimension > 1)
    errY = H5PartReadDataFloat64(dataFile,"y",dataY);
  if(gdimension > 2)
    errZ = H5PartReadDataFloat64(dataFile,"z",dataZ);

  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read x" << std::endl;
  if(errY != H5_SUCCESS) 
    rank || clog(warn) << "Unable to read y" << std::endl;
  if(errZ != H5_SUCCESS) 
    rank || clog(warn) << "Unable to read z" << std::endl;


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
  errX = H5PartReadDataFloat64(dataFile,"vx",dataX);
  if(gdimension > 1)
    errY = H5PartReadDataFloat64(dataFile,"vy",dataY);
  if(gdimension > 2)
    errZ = H5PartReadDataFloat64(dataFile,"vz",dataZ);
 
  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read vx" << std::endl;
  if(errY != H5_SUCCESS) 
    rank || clog(warn) << "Unable to read vy" << std::endl;
  if(errZ != H5_SUCCESS) 
    rank || clog(warn) << "Unable to read vz" << std::endl;

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
  errX = H5PartReadDataFloat64(dataFile,"ax",dataX);
  if(gdimension > 1)
    errY = H5PartReadDataFloat64(dataFile,"ay",dataY);
  if(gdimension > 2)
    errZ = H5PartReadDataFloat64(dataFile,"az",dataZ);
  
  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read ax" << std::endl;
  if(errY != H5_SUCCESS) 
    rank || clog(warn) << "Unable to read ay" << std::endl;
  if(errZ != H5_SUCCESS) 
    rank || clog(warn) << "Unable to read az" << std::endl;

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

  errX = H5PartReadDataFloat64(dataFile,"m",dataX); 
  errY = H5PartReadDataFloat64(dataFile,"rho",dataY);
  errZ = H5PartReadDataFloat64(dataFile,"h",dataZ);

  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read m" << std::endl;
  if(errY != H5_SUCCESS)
    rank || clog(warn) << "Unable to read rho" << std::endl;
  if(errZ != H5_SUCCESS)
    rank || clog(warn) << "Unable to read h" << std::endl;

  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setMass(dataX[i]);
    bodies[i].second.setDensity(dataY[i]);
    bodies[i].second.setSmoothinglength(dataZ[i]);
  }

  // Reset buffer to 0, if next value not present 
  std::fill(dataX,dataX+nparticlesproc,0.);
  std::fill(dataY,dataY+nparticlesproc,0.);
  std::fill(dataZ,dataZ+nparticlesproc,0.); 

  // Pressure  
  errX = H5PartReadDataFloat64(dataFile,"P",dataX);
  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read P" <<std::endl;
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setPressure(dataX[i]);
  }

  // Internal Energy  
  #ifdef INTERNAL_ENERGY
  rank|| clog(trace)<<"Reading internal energy"<<std::endl;
  std::fill(dataX,dataX+nparticlesproc,0.);
  errX = H5PartReadDataFloat64(dataFile,"u",dataX);
  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read u"<<std::endl;
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
    rank|| clog(trace)<<"Setting ID for particles"<<std::endl;
    // Otherwise generate the id 
    int64_t start = (totalnbodies/size)*rank+1;
    for(int64_t i=0; i<nparticlesproc; ++i){
      bodies[i].second.setId(start+i); 
    }
    if(size == rank + 1){
      assert(totalnbodies==bodies.back().second.getId()); 
    }
  }
  
  // Reset buffer to 0, if next value not present 
  std::fill(dataX,dataX+nparticlesproc,0.);

  // delta T   
  errX = H5PartReadDataFloat64(dataFile,"dt",dataX);
  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read dt" << std::endl;
  for(int64_t i=0; i<nparticlesproc; ++i){
    bodies[i].second.setDt(dataX[i]);
  }
  
  // Reset buffer to 0, if next value not present 
  std::fill(dataInt32,dataInt32+nparticlesproc,0.);

  // delta T   
  errX = H5PartReadDataInt32(dataFile,"type",dataInt32);
  if(errX != H5_SUCCESS)
    rank || clog(warn) << "Unable to read type"<< std::endl;
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

  rank|| clog(trace)<<"Input particles.done"<<std::endl;

  // Restore the H5HUT error level
  H5SetVerbosityLevel  (1); 
}// inputDataHDF5

// Output data in HDF5 format 
// Generate the associate XDMF file 
void outputDataHDF5(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    const char* fileprefix,
    int64_t iteration,
    double totaltime)
{
  int step = iteration/param::out_h5data_every; 

  bool do_diff_files = param::out_h5data_separate_iterations;
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  MPI_Barrier(MPI_COMM_WORLD);
  rank|| clog(trace)<<"Output particles"<<std::flush;

  int64_t nparticlesproc = bodies.size();

  char filename[128];
  if(do_diff_files){
    sprintf(filename,"%s_%05d.h5part",fileprefix,step);
  }
  else {
    sprintf(filename,"%s.h5part",fileprefix);
  }

  // If different file per iteration, just remove the file with same name 
  // If one big file, remove the file if it is the Step 0  
  if(do_diff_files && rank == 0){
    remove(filename);
  }else if(step == 0 && rank == 0){
    remove(filename);
  }

  // Wait for removing the file before writing in 
  MPI_Barrier(MPI_COMM_WORLD);

  h5_file_t* dataFile = H5OpenFile(filename,H5_O_RDWR | H5_VFD_MPIIO_IND,
      MPI_COMM_WORLD);
  
  //-------------------GLOBAL HEADER-------------------------------------------
  // Only for the first output
  if(do_diff_files || step == 0){
    int gdimension32 = gdimension;  
    H5WriteFileAttribInt32(dataFile,"dimension",&gdimension32,1);
  }

  //------------------STEP HEADER----------------------------------------------
  // Put the step header
  H5SetStep(dataFile,step);
  H5WriteStepAttribFloat64(dataFile,"time",&physics::totaltime,1); 
  H5WriteStepAttribInt64(dataFile,"iteration",&physics::iteration,1);
  H5WriteStepAttribFloat64(dataFile,"timestep",&physics::dt,1);  
  //------------------STEP DATA------------------------------------------------

  H5PartSetNumParticles(dataFile,nparticlesproc);

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
  H5PartWriteDataFloat64(dataFile,"x",b1); 
  H5PartWriteDataFloat64(dataFile,"y",b2);
  H5PartWriteDataFloat64(dataFile,"z",b3);

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
  H5PartWriteDataFloat64(dataFile,"vx",b1);
  H5PartWriteDataFloat64(dataFile,"vy",b2);
  H5PartWriteDataFloat64(dataFile,"vz",b3);

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
  H5PartWriteDataFloat64(dataFile,"ax",b1);
  H5PartWriteDataFloat64(dataFile,"ay",b2);
  H5PartWriteDataFloat64(dataFile,"az",b3);

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
  H5PartWriteDataFloat64(dataFile,"h",b1);
  H5PartWriteDataFloat64(dataFile,"rho",b2);
#ifdef INTERNAL_ENERGY
  H5PartWriteDataFloat64(dataFile,"u",b3);
#endif 

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
  H5PartWriteDataFloat64(dataFile,"P",b1);
  H5PartWriteDataFloat64(dataFile,"m",b2);
  H5PartWriteDataFloat64(dataFile,"dt",b3);
  H5PartWriteDataInt64(dataFile,"id",bi);
  H5PartWriteDataInt32(dataFile,"type",bint);

  // Output the rank for analysis 
  std::fill(bi,bi+nparticlesproc,rank);
  H5PartWriteDataInt64(dataFile,"rank",bi);

  H5CloseFile(dataFile);

  delete[] b1;
  delete[] b2;
  delete[] b3;
  delete[] bi;
  delete[] bint;

  rank|| clog(trace)<<".done"<<std::endl;
  
}// outputDataHDF5

} // namespace io

#endif


