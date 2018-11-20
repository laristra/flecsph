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

#include <hdf5.h>
//#include <H5hut.h>

//#include "tree.h"
//#include "physics.h"

namespace io{

hid_t IO_group_id; // Group id to keep track of the current step
// Data for hyperslab
hsize_t IO_offset;
hsize_t IO_count;

int64_t IO_nparticlesproc;
int64_t IO_nparticles;

hid_t
H5P_openFile(const char * filename,	unsigned int flags )
{
  //std::cout<<"Opening file";
  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;
  /* Set up file access property list with parallel I/O access */
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);

  hid_t file_id = 0;
  if( access( filename, F_OK ) != -1 ) {
    // file exists
    //std::cout<<"EXISTING FILE"<<std::endl;
    file_id = H5Fopen(filename, flags, plist_id);
  }else{
    //std::cout<<"CREATING FILE"<<std::endl;
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  }

  H5Pclose(plist_id);
  //std::cout<<"done"<<std::endl;
  return file_id;
}

void
H5P_closeFile(hid_t& file_id)
{
  //if(H5I_type_t(group_id) == H5I_GROUP)
  H5Gclose(IO_group_id);
  H5Fclose(file_id);
}

bool
H5P_hasStep(hid_t& file_id, size_t step)
{
  hid_t error_stack = 0;
  H5Eset_auto(error_stack, NULL, NULL);

  /* Turn off error handling */
  //H5Eset_auto(error_stack, NULL, NULL);

  /* Save old error handler */
  //herr_t (*old_func)(hid_t,void*);
  //void *old_client_data;
  //hid_t estack_id = H5Eget_current_stack();
  //hid_t nstack = H5Ecreate_stack();
  //H5Eget_auto(estack_id, &old_func, &old_client_data);
  /* Turn off error handling */
  //H5Eset_auto(nstack, NULL, NULL);
  //std::cout<<"Checking step: "<<step<<" ";
  char cstep[255];
  sprintf(cstep,"/Step#%lu",step);
  hid_t stat = H5Gget_objinfo (file_id, cstep, 0, NULL);
  /* Restore previous error handler */
  //H5Eset_auto(estack_id, old_func, old_client_data);
  if (stat == 0){
    //std::cout<<"FOUND"<<std::endl;
    return true;
  }
  //std::cout<<"NOTFOUND"<<std::endl;
  return false;
}

void
H5P_setStep(
  hid_t& file_id,
  size_t step
)
{
  //std::cout<<"Creating GROUP";
  char cstep[255];
  sprintf(cstep,"/Step#%lu",step);
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  if(H5P_hasStep(file_id,step)){
    IO_group_id = H5Gopen(file_id, cstep, H5P_DEFAULT );
  }else{
    IO_group_id = H5Gcreate(file_id, cstep, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  //std::cout<<"done"<<std::endl;
  H5Pclose(plist_id);
}

/**
*
*/
template<
  typename T>
hid_t
H5P_writeDataset(
  hid_t& file_id,
  const char * dsname,
  T* data,
  size_t dim = IO_nparticlesproc)
{

  hid_t type = H5T_NATIVE_INT;
  if(typeid(T) == typeid(int)){
  }else if(typeid(T) == typeid(double)){
    type = H5T_NATIVE_DOUBLE;
  }else if (typeid(T) == typeid(int64_t)){
    type = H5T_NATIVE_LLONG;
  }else if (typeid(T) == typeid(uint64_t)){
    type = H5T_NATIVE_ULLONG;
  } else {
    std::cout<<"Unknown type: "<<typeid(T).name()<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  //std::cout<<"Writing Dataset "<<dsname<<std::flush;
  hid_t status = 1;
  //char cstep[255];
  //sprintf(cstep,"/#Step%lu/%s",step,dsname);
  hsize_t hdim = dim;
  /* Create the dataspace for the dataset.*/
  hsize_t total = IO_nparticles;
  //std::cout<<"Total: "<<total<<std::endl;
  hid_t filespace = H5Screate_simple(1, &total, NULL);
  hid_t dset_id = H5Dcreate(IO_group_id, dsname, type, filespace,
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);

  hsize_t offset_in = 0;
  hsize_t count_in = dim;
  //std::cout<<"Count: "<<dim<<std::endl;
  hid_t memspace = H5Screate_simple(1, &count_in, NULL);
  // output hyperslab
  //status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &offset_in, NULL,
  //   &count_in, NULL);

  // Select the hyperslab
  hid_t dataspace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &IO_offset, NULL, &IO_count, NULL);

  /*Create property list for collective dataset write.*/
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  status = H5Dwrite(dset_id, type, memspace, dataspace, plist_id, data);

  // Close everythings
  H5Sclose(memspace);
  H5Sclose(dataspace);
  //H5Sclose(filespace);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);

  //std::cout<<"done"<<std::endl;
  return status;
}


template<
  typename T>
hid_t
H5P_readAttribute(
  hid_t& file_id,
  const char * dsname,
  T* data)
{
  //std::cout<<"Reading Attribute";
  //char cstep[255];
  //sprintf(cstep,"/%s",dsname);
  hid_t status = 1;
  hid_t att_id = H5Aopen(file_id,dsname,H5P_DEFAULT);
  /*Create property list for collective dataset write.*/
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  /* To write dataset independently use
  * H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);*/
  if(typeid(T) == typeid(int)){
    status = H5Aread(att_id, H5T_NATIVE_INT, data);
  }else if(typeid(T) == typeid(double)){
    status = H5Aread(att_id, H5T_NATIVE_DOUBLE, data);
  }else if (typeid(T) == typeid(int64_t)){
    status = H5Aread(att_id, H5T_NATIVE_LLONG, data);
  }else if (typeid(T) == typeid(uint64_t)){
    status = H5Aread(att_id, H5T_NATIVE_ULLONG, data);
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  H5Aclose(att_id);
  H5Pclose(plist_id);
  //std::cout<<"done: "<<*data<<std::endl;
  return status;
}

template<
  typename T>
hid_t
H5P_readAttributeStep(
  hid_t& file_id,
  const char * dsname,
  T* data)
{
  //std::cout<<"Reading Attribute Step: "<<dsname;
  hid_t status = 1;
  hid_t att_id = H5Aopen(IO_group_id,dsname,H5P_DEFAULT);
  /*Create property list for collective dataset write.*/
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  /* To write dataset independently use
  * H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);*/
  if(typeid(T) == typeid(int)){
    status = H5Aread(att_id, H5T_NATIVE_INT, data);
  }else if(typeid(T) == typeid(double)){
    status = H5Aread(att_id, H5T_NATIVE_DOUBLE, data);
  }else if (typeid(T) == typeid(int64_t)){
    status = H5Aread(att_id, H5T_NATIVE_LLONG, data);
  }else if (typeid(T) == typeid(uint64_t)){
    status = H5Aread(att_id, H5T_NATIVE_ULLONG, data);
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  H5Aclose(att_id);
  H5Pclose(plist_id);
  //std::cout<<"done: "<<*data<<std::endl;
  return status;
}

template<
  typename T>
hid_t
H5P_writeAttribute(
  hid_t& file_id,
  const char * dsname,
  T* data)
{
  //std::cout<<"Writing Attribute";
  hid_t status = 1;

  hsize_t hdim = 1;
  /* Create the dataspace for the dataset.*/
  hid_t filespace = H5Screate_simple(1, &hdim, NULL);
  /*Create the dataset with default properties and close filespace.*/
  hid_t att_id = 1;

  if(typeid(T) == typeid(int)){
    att_id = H5Acreate(file_id, dsname, H5T_NATIVE_INT, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_INT, data);
  }else if(typeid(T) == typeid(double)){
    att_id = H5Acreate(file_id, dsname, H5T_NATIVE_DOUBLE, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_DOUBLE, data);
  }else if (typeid(T) == typeid(int64_t)){
    att_id = H5Acreate(file_id, dsname, H5T_NATIVE_LLONG, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_LLONG, data);
  }else if (typeid(T) == typeid(uint64_t)){
    att_id = H5Acreate(file_id, dsname, H5T_NATIVE_ULLONG, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_ULLONG, data);
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  H5Sclose(filespace);
  H5Aclose(att_id);
  //std::cout<<"done"<<std::endl;
  return status;
}

template<
  typename T>
hid_t
H5P_writeAttributeStep(
  hid_t& file_id,
  const char * dsname,
  T* data)
{
  //std::cout<<"Writing Attribute Step";
  hid_t status = 1;
  hsize_t hdim = 1;
  //char cstep[255];
  //sprintf(cstep,"/#Step%lu/%s",step,dsname);
  /* Create the dataspace for the dataset.*/
  hid_t filespace = H5Screate_simple(1, &hdim, NULL);
  /*Create the dataset with default properties and close filespace.*/
  hid_t att_id = 1;

  if(typeid(T) == typeid(int)){
    att_id = H5Acreate(IO_group_id, dsname, H5T_NATIVE_INT, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_INT, data);
  }else if(typeid(T) == typeid(double)){
    att_id = H5Acreate(IO_group_id, dsname, H5T_NATIVE_DOUBLE, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_DOUBLE, data);
  }else if (typeid(T) == typeid(int64_t)){
    att_id = H5Acreate(IO_group_id, dsname, H5T_NATIVE_LLONG, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_LLONG, data);
  }else if (typeid(T) == typeid(uint64_t)){
    att_id = H5Acreate(IO_group_id, dsname, H5T_NATIVE_ULLONG, filespace,
      H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(att_id, H5T_NATIVE_ULLONG, data);
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  H5Sclose(filespace);
  H5Aclose(att_id);
  //std::cout<<"done"<<std::endl;
  return status;
}

template<
  typename T>
hid_t
H5P_readDataset(
  hid_t& file_id,
  const char * dsname,
  T* data,
  size_t dim = IO_nparticlesproc)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //std::cout<<"Reading Dataset: "<<dsname<<std::flush;
  hid_t status = 1;
  /* Open the dataset*/
  hid_t dset_id = H5Dopen (IO_group_id, dsname, H5P_DEFAULT);
  /*Create property list for collective dataset write.*/
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  hsize_t out = dim;
  hsize_t offset_out = 0;
  hsize_t count_out = dim;
  hid_t memspace = H5Screate_simple(1, &out, NULL);
  // output hyperslab
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &offset_out, NULL,
     &count_out, NULL);

  // Select the hyperslab
  hid_t dataspace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &IO_offset, NULL, &IO_count, NULL);
  /* To write dataset independently use
  * H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);*/
  if(typeid(T) == typeid(int)){
      status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, dataspace,
        plist_id, data);
  }else if(typeid(T) == typeid(double)){
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, dataspace,
        plist_id, data);
  }else if (typeid(T) == typeid(int64_t)){
    status = H5Dread(dset_id, H5T_NATIVE_LLONG, memspace, dataspace,
      plist_id, data);
  }else if (typeid(T) == typeid(uint64_t)){
    status = H5Dread(dset_id, H5T_NATIVE_ULLONG, memspace, dataspace,
      plist_id, data);
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  //std::cout<<"done"<<std::endl<<std::flush;
  return status;
}

size_t
H5P_setNumParticles(const int64_t& nparticlesproc)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  // Compute arrays for hyperslab

  IO_nparticlesproc = nparticlesproc;
  std::vector<int64_t> offcount(size+1);
  MPI_Allgather(&nparticlesproc,1,MPI_INT64_T,
    &offcount[0],1,MPI_INT64_T,MPI_COMM_WORLD);

  for(int i = 1 ; i < size; ++i)
    offcount[i] += offcount[i-1];
  offcount.insert(offcount.begin(),0);

  IO_nparticles = offcount[size];
  IO_offset = offcount[rank];
  IO_count = nparticlesproc;

  //std::cout<<"nparticles: "<<IO_nparticles<<" nparticlesproc: "<<IO_nparticlesproc<<std::endl;
  return IO_nparticlesproc;
}

size_t
H5P_getNumParticles(hid_t file_id)
{
  //std::cout<<"Get N PART";
  // Read x and get the dimension of the dataset
  hid_t status = 1 ;
  //char cstep[255];
  //sprintf(cstep,"/#Step%lu/x",step);
  hid_t dset_id = H5Dopen (IO_group_id, "x", H5P_DEFAULT);
  hid_t dspace = H5Dget_space(dset_id);
  //size_t parts = H5Sget_simple_extent_ndims(dspace);
  hsize_t dim;
  H5Sget_simple_extent_dims(dspace, &dim, NULL);
  size_t parts = dim;
  H5Dclose(dset_id);
  H5Sclose(dspace);
  //std::cout<<"done: "<<parts<<std::endl;
  return parts;
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
  //H5SetVerbosityLevel  (0);

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  rank|| clog(trace)<<"Input particles" << std::endl;

  hid_t dataFile = H5P_openFile(filename,H5F_ACC_RDONLY);

  // The input file depend of the type of output, separate or one file.
  //h5_file_t * dataFile = nullptr;
  //if(startIteration == 0){
  //  dataFile = H5OpenFile(filename,H5_O_RDONLY
  //    | H5_VFD_MPIIO_IND,// Flag to be not use mpiposix
  //    MPI_COMM_WORLD);
  //}

  // ------------- CHECK IF THE SEARCHED ITERATION EXISTS ---------------------
  // Go through all the steps of the file and try to read the iteration
  if(!param::out_h5data_separate_iterations && startIteration != 0){
   H5P_closeFile(dataFile);
   dataFile = H5P_openFile(filename,H5F_ACC_RDONLY);
   int step = 1;
    bool end = false;
    bool found = false;
    while(!end){
      int hasStep = H5P_hasStep(dataFile,step);
      if(hasStep){
        // Check iteration number
        H5P_setStep(dataFile,step);
        int64_t iteration;
        if(0 != H5P_readAttributeStep(dataFile,"iteration",&iteration)){
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
      auto stepFile = H5P_openFile(step_filename,H5F_ACC_RDONLY);
      int64_t iteration;
      H5P_setStep(dataFile,step);
      if(0 != H5P_readAttributeStep(stepFile,"iteration",&iteration)){
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
      H5P_closeFile(stepFile);
    }
    startStep = step-1;
    char diff_filename[128];
    sprintf(diff_filename,"%s_%05d.h5part",output_file_prefix,startStep);
    rank || clog(warn) << "Reading from file "<<diff_filename<<std::endl;
    // Change input file in this case
    dataFile = H5P_openFile(diff_filename,H5F_ACC_RDONLY);
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
      auto outputFile = H5P_openFile(output_filename,H5F_ACC_RDONLY);
      H5P_setStep(outputFile,startStep);
      // Check if startStep == lastStep
      int lastStep = startStep;
      while(H5P_hasStep(outputFile,lastStep)){lastStep++;}
      --lastStep;
      rank || clog(trace)<<"startStep: "<<startStep<<" lastStep: "
        <<lastStep<<std::endl;
      if(startStep != lastStep){
        rank || clog(error) << "First step not last step in output"<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
      }
      H5P_closeFile(outputFile);
    }
  }


  if(dataFile == 0){
    rank || clog(error) << "Cannot find data file"<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }

  int val = H5P_hasStep(dataFile,startStep);

  // If yes, get into this step
  H5P_setStep(dataFile,startStep);
  // Get the number of particles
  int64_t nparticles = 0UL;
  nparticles = H5P_getNumParticles(dataFile);

  int64_t nparticlesproc = nparticles/size;
  // Handle the number of particles for the last one
  if(size == rank + 1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  H5P_setNumParticles(nparticlesproc);

  // Register for main
  totalnbodies = nparticles;
  nbodies = IO_nparticlesproc;

  //--------------- READ GLOBAL ATTRIBUTES ------------------------------------
  // read the number of dimension
  int32_t dimension;
  if(0 == H5P_readAttribute(dataFile,"dimension",&dimension)){
    assert(gdimension == dimension);
  }else{
    rank|| clog(error)<<"No dimension value: setting to default 3"<<std::endl;
    dimension = 3;
  }

  //--------------- READ DATA FROM STEP ---------------------------------------
  // Set the number of particles read by each process
  //H5PartSetNumParticles(dataFile,nparticlesproc);
  // Resize the bodies vector
  bodies.clear();
  bodies.resize(IO_nparticlesproc);

  // Try to read timestep if exists
  if(startIteration != 0){
    double timestep = 0.;
    double totaltime = 0.;
    if(0 == H5P_readAttributeStep(dataFile,"timestep",&timestep)){
      physics::dt = timestep;
    }else{
      rank || clog(warn)<<"Unable to read timestep from file"<<std::endl;
    }
    if(0 == H5P_readAttributeStep(dataFile,"time",&totaltime)){
      physics::totaltime = totaltime;
    }else{
      rank || clog(warn)<<"Unable to read totaltime from file"<<std::endl;
    }
  }

  // Read the dataset and fill the particles data
  double* dataX = new double[IO_nparticlesproc];
  double* dataY = new double[IO_nparticlesproc];
  double* dataZ = new double[IO_nparticlesproc];
  int64_t* dataInt = new int64_t[IO_nparticlesproc];
  int * dataInt32 = new int[IO_nparticlesproc];

  // Handle errors from H5HUT
  hid_t errX, errY, errZ;
  errX = errY = errZ = 0; // prevent warnings

  // Positions
  errX = H5P_readDataset(dataFile,"x",dataX);
  if(gdimension > 1)
    errY = H5P_readDataset(dataFile,"y",dataY);
  if(gdimension > 2)
    errZ = H5P_readDataset(dataFile,"z",dataZ);

  if(errX != 0)
    rank || clog(warn) << "Unable to read x" << std::endl;
  if(errY != 0)
    rank || clog(warn) << "Unable to read y" << std::endl;
  if(errZ != 0)
    rank || clog(warn) << "Unable to read z" << std::endl;


  for(int64_t i=0; i<IO_nparticlesproc; ++i){
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
  std::fill(dataX,dataX+IO_nparticlesproc,0.);
  std::fill(dataY,dataY+IO_nparticlesproc,0.);
  std::fill(dataZ,dataZ+IO_nparticlesproc,0.);

  // Velocity
  errX = H5P_readDataset(dataFile,"vx",dataX);
  if(gdimension > 1)
    errY = H5P_readDataset(dataFile,"vy",dataY);
  if(gdimension > 2)
    errZ = H5P_readDataset(dataFile,"vz",dataZ);

  if(errX != 0)
    rank || clog(warn) << "Unable to read vx" << std::endl;
  if(errY != 0)
    rank || clog(warn) << "Unable to read vy" << std::endl;
  if(errZ != 0)
    rank || clog(warn) << "Unable to read vz" << std::endl;

  for(int64_t i=0; i<IO_nparticlesproc; ++i){
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
  std::fill(dataX,dataX+IO_nparticlesproc,0.);
  std::fill(dataY,dataY+IO_nparticlesproc,0.);
  std::fill(dataZ,dataZ+IO_nparticlesproc,0.);

  // Acceleration
  errX = H5P_readDataset(dataFile,"ax",dataX);
  if(gdimension > 1)
    errY = H5P_readDataset(dataFile,"ay",dataY);
  if(gdimension > 2)
    errZ = H5P_readDataset(dataFile,"az",dataZ);

  if(errX != 0)
    rank || clog(warn) << "Unable to read ax" << std::endl;
  if(errY != 0)
    rank || clog(warn) << "Unable to read ay" << std::endl;
  if(errZ != 0)
    rank || clog(warn) << "Unable to read az" << std::endl;

  for(int64_t i=0; i<IO_nparticlesproc; ++i){
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
  std::fill(dataX,dataX+IO_nparticlesproc,0.);
  std::fill(dataY,dataY+IO_nparticlesproc,0.);
  std::fill(dataZ,dataZ+IO_nparticlesproc,0.);

  errX = H5P_readDataset(dataFile,"m",dataX);
  errY = H5P_readDataset(dataFile,"rho",dataY);
  errZ = H5P_readDataset(dataFile,"h",dataZ);

  if(errX != 0)
    rank || clog(warn) << "Unable to read m" << std::endl;
  if(errY != 0)
    rank || clog(warn) << "Unable to read rho" << std::endl;
  if(errZ != 0)
    rank || clog(warn) << "Unable to read h" << std::endl;

  for(int64_t i=0; i<IO_nparticlesproc; ++i){
    bodies[i].second.setMass(dataX[i]);
    bodies[i].second.setDensity(dataY[i]);
    bodies[i].second.setSmoothinglength(dataZ[i]);
  }

  // Reset buffer to 0, if next value not present
  std::fill(dataX,dataX+IO_nparticlesproc,0.);
  std::fill(dataY,dataY+IO_nparticlesproc,0.);
  std::fill(dataZ,dataZ+IO_nparticlesproc,0.);

  // Pressure
  errX = H5P_readDataset(dataFile,"P",dataX);
  if(errX != 0)
    rank || clog(warn) << "Unable to read P" <<std::endl;
  for(int64_t i=0; i<IO_nparticlesproc; ++i){
    bodies[i].second.setPressure(dataX[i]);
  }

  // Internal Energy
  #ifdef INTERNAL_ENERGY
  rank|| clog(trace)<<"Reading internal energy"<<std::endl;
  std::fill(dataX,dataX+IO_nparticlesproc,0.);
  errX = H5P_readDataset(dataFile,"u",dataX);
  if(errX != 0)
    rank || clog(warn) << "Unable to read u"<<std::endl;
  for(int64_t i=0; i<IO_nparticlesproc; ++i){
    bodies[i].second.setInternalenergy(dataX[i]);
  }
  #endif



  //bool b_index = false;
  // \TODO check if user ID is uniq
  bool b_index = 0 ==
  H5P_readDataset(dataFile,"id",dataInt);
  if(b_index){
    for(int64_t i=0; i<IO_nparticlesproc; ++i){
      bodies[i].second.setId(dataInt[i]);
    }
  }else{
    rank|| clog(trace)<<"Setting ID for particles"<<std::endl;
    // Otherwise generate the id
    int64_t start = (totalnbodies/size)*rank+1;
    for(int64_t i=0; i<IO_nparticlesproc; ++i){
      bodies[i].second.setId(start+i);
    }
    if(size == rank + 1){
      assert(totalnbodies==bodies.back().second.getId());
    }
  }

  // Reset buffer to 0, if next value not present
  std::fill(dataX,dataX+IO_nparticlesproc,0.);

  // delta T
  errX = H5P_readDataset(dataFile,"dt",dataX);
  if(errX != 0)
    rank || clog(warn) << "Unable to read dt" << std::endl;
  for(int64_t i=0; i<IO_nparticlesproc; ++i){
    bodies[i].second.setDt(dataX[i]);
  }

  // Reset buffer to 0, if next value not present
  std::fill(dataInt32,dataInt32+IO_nparticlesproc,0.);

  // delta T
  errX = H5P_readDataset(dataFile,"type",dataInt32);
  if(errX != 0)
    rank || clog(warn) << "Unable to read type"<< std::endl;
  for(int64_t i=0; i<IO_nparticlesproc; ++i){
    bodies[i].second.setType(dataInt32[i]);
  }

  delete[] dataX;
  delete[] dataY;
  delete[] dataZ;
  delete[] dataInt;
  delete[] dataInt32;

  MPI_Barrier(MPI_COMM_WORLD);
  H5Fclose(dataFile);
  //H5CloseFile(dataFile);

  rank|| clog(trace)<<"Input particles.done"<<std::endl;


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
  // Check if file exists
  hid_t dataFile = H5P_openFile(filename,H5F_ACC_RDWR);

  //-------------------GLOBAL HEADER-------------------------------------------
  // Only for the first output
  if(do_diff_files || step == 0){
    int gdimension32 = gdimension;
    H5P_writeAttribute(dataFile,"dimension",&gdimension32);
  }

  //------------------STEP HEADER----------------------------------------------
  // Put the step header
  H5P_setStep(dataFile,step);
  H5P_writeAttributeStep(dataFile,"time",&physics::totaltime);
  H5P_writeAttributeStep(dataFile,"iteration",&physics::iteration);
  H5P_writeAttributeStep(dataFile,"timestep",&physics::dt);
  //------------------STEP DATA------------------------------------------------

  H5P_setNumParticles(bodies.size());

  // 4 buffers, 3 double and 1 int64
  double* b1 = new double[IO_nparticlesproc];
  double* b2 = new double[IO_nparticlesproc];
  double* b3 = new double[IO_nparticlesproc];
  int64_t* bi = new int64_t[IO_nparticlesproc];
  int32_t* bint = new int32_t[IO_nparticlesproc];

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
  H5P_writeDataset(dataFile,"x",b1);
  H5P_writeDataset(dataFile,"y",b2);
  H5P_writeDataset(dataFile,"z",b3);


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
  H5P_writeDataset(dataFile,"vx",b1);
  H5P_writeDataset(dataFile,"vy",b2);
  H5P_writeDataset(dataFile,"vz",b3);

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
  H5P_writeDataset(dataFile,"ax",b1);
  H5P_writeDataset(dataFile,"ay",b2);
  H5P_writeDataset(dataFile,"az",b3);

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
  H5P_writeDataset(dataFile,"h",b1);
  H5P_writeDataset(dataFile,"rho",b2);
#ifdef INTERNAL_ENERGY
  H5P_writeDataset(dataFile,"u",b3);
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
  H5P_writeDataset(dataFile,"P",b1);
  H5P_writeDataset(dataFile,"m",b2);
  H5P_writeDataset(dataFile,"dt",b3);
  H5P_writeDataset(dataFile,"id",bi);
  H5P_writeDataset(dataFile,"type",bint);

  // Output the rank for analysis
  std::fill(bi,bi+IO_nparticlesproc,rank);
  H5P_writeDataset(dataFile,"rank",bi);

  H5P_closeFile(dataFile);

  delete[] b1;
  delete[] b2;
  delete[] b3;
  delete[] bi;
  delete[] bint;

  rank|| clog(trace)<<".done"<<std::endl;

}// outputDataHDF5
} // namespace io

#endif // _mpisph_io_h_
