/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
#include <dirent.h>
#include <iostream>
#include <libgen.h>
#include <vector>

#include "default_physics.h"
#include "params.h"

#include <hdf5.h>
//#include <H5hut.h>

//#include "tree.h"
//#include "physics.h"

namespace io {

hid_t IO_group_id; // Group id to keep track of the current step
// Data for hyperslab
hsize_t IO_offset;
hsize_t IO_count;
static int output_step = 0;
const int MAX_FNAME_LEN = 256;
// TODO: overload ostream instead, i.e.smth like, clog_exit << "ERROR!"
#define FULLSTOP exit(MPI_Barrier(MPI_COMM_WORLD) && MPI_Finalize());

int64_t IO_nparticlesproc;
int64_t IO_nparticles;

MPI_Comm comm_ = MPI_COMM_WORLD;

template <typename T> hid_t H5P_getType(T *data) {
  hid_t type = H5T_NATIVE_INT;
  if (typeid(T) == typeid(int)) {
  } else if (typeid(T) == typeid(double)) {
    type = H5T_NATIVE_DOUBLE;
  } else if (typeid(T) == typeid(int64_t)) {
    type = H5T_NATIVE_LLONG;
  } else if (typeid(T) == typeid(uint64_t)) {
    type = H5T_NATIVE_ULLONG;
  } else {
    std::cout << "Unknown type: " << typeid(T).name() << std::endl;
    MPI_Barrier(comm_);
    MPI_Finalize();
  }
  return type;
}

hid_t H5P_openFile(const char *filename, unsigned int flags) {
  MPI_Comm comm = comm_;
  MPI_Info info = MPI_INFO_NULL;
  /* Set up file access property list with parallel I/O access */
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);

  hid_t file_id = 0;
  if (access(filename, F_OK) != -1) {
    file_id = H5Fopen(filename, flags, plist_id);
  } else {
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  }

  H5Pclose(plist_id);
  return file_id;
}

void H5P_closeFile(hid_t &file_id) {
  H5Gclose(IO_group_id);
  H5Fclose(file_id);
}

bool H5P_hasStep(hid_t &file_id, size_t step) {
  hid_t error_stack = 0;
  H5Eset_auto(error_stack, NULL, NULL);

  char cstep[255];
  sprintf(cstep, "/Step#%lu", step);
  hid_t stat = H5Gget_objinfo(file_id, cstep, 0, NULL);
  return !(stat); // true if found (stat==0), false if not
}

void H5P_setStep(hid_t &file_id, size_t step) {
  char cstep[255];
  sprintf(cstep, "/Step#%lu", step);
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  if (H5P_hasStep(file_id, step)) {
    IO_group_id = H5Gopen(file_id, cstep, H5P_DEFAULT);
  } else {
    IO_group_id =
        H5Gcreate(file_id, cstep, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  H5Pclose(plist_id);
}

template <typename T>
hid_t H5P_writeDataset(hid_t &file_id, const char *dsname, T *data,
                       size_t dim = IO_nparticlesproc) {

  hid_t type = H5P_getType(data);

  hid_t status = 1;
  hsize_t hdim = dim;
  /* Create the dataspace for the dataset.*/
  hsize_t total = IO_nparticles;
  hid_t filespace = H5Screate_simple(1, &total, NULL);
  hid_t dset_id = H5Dcreate(IO_group_id, dsname, type, filespace, H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);

  hsize_t offset_in = 0;
  hsize_t count_in = dim;
  hid_t memspace = H5Screate_simple(1, &count_in, NULL);

  // Select the hyperslab
  hid_t dataspace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &IO_offset, NULL, &IO_count,
                      NULL);

  /*Create property list for collective dataset write.*/
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  status = H5Dwrite(dset_id, type, memspace, dataspace, plist_id, data);

  // Close everythings
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  return status;
}

template <typename T>
hid_t H5P_readAttribute(hid_t &file_id, const char *dsname, T *data) {
  hid_t type = H5P_getType(data);

  hid_t status = 1;
  hid_t att_id = H5Aopen(file_id, dsname, H5P_DEFAULT);
  /*Create property list for collective dataset write.*/
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  status = H5Aread(att_id, type, data);

  H5Aclose(att_id);
  H5Pclose(plist_id);
  return status;
}

template <typename T>
hid_t H5P_readAttributeStep(hid_t &file_id, const char *dsname, T *data) {
  hid_t type = H5P_getType(data);

  hid_t status = 1;
  hid_t att_id = H5Aopen(IO_group_id, dsname, H5P_DEFAULT);
  /*Create property list for collective dataset write.*/
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  status = H5Aread(att_id, type, data);

  H5Aclose(att_id);
  H5Pclose(plist_id);
  return status;
}

template <typename T>
hid_t H5P_writeAttribute(hid_t &file_id, const char *dsname, T *data) {
  hid_t type = H5P_getType(data);

  hid_t status = 1;
  hsize_t hdim = 1;
  /* Create the dataspace for the dataset.*/
  hid_t filespace = H5Screate_simple(1, &hdim, NULL);
  /*Create the dataset with default properties and close filespace.*/
  hid_t att_id = 1;

  att_id =
      H5Acreate(file_id, dsname, type, filespace, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att_id, type, data);

  H5Sclose(filespace);
  H5Aclose(att_id);
  return status;
}

template <typename T>
hid_t H5P_writeAttributeStep(hid_t &file_id, const char *dsname, T *data) {
  hid_t type = H5P_getType(data);

  hid_t status = 1;
  hsize_t hdim = 1;
  hid_t filespace = H5Screate_simple(1, &hdim, NULL);
  /*Create the dataset with default properties and close filespace.*/
  hid_t att_id = 1;

  att_id =
      H5Acreate(IO_group_id, dsname, type, filespace, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att_id, type, data);

  H5Sclose(filespace);
  H5Aclose(att_id);
  return status;
}

template <typename T>
hid_t H5P_readDataset(hid_t &file_id, const char *dsname, T *data,
                      size_t dim = IO_nparticlesproc) {
  hid_t type = H5P_getType(data);

  int rank;
  MPI_Comm_rank(comm_, &rank);
  hid_t status = 1;
  /* Open the dataset*/
  hid_t dset_id = H5Dopen(IO_group_id, dsname, H5P_DEFAULT);
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
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &IO_offset, NULL, &IO_count,
                      NULL);

  status = H5Dread(dset_id, type, memspace, dataspace, plist_id, data);

  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  return status;
}

/*
 * Read scalar data from HDF5 file and assign to bodies
 *  - if dataset doesn't exist, produce warning;
 *  - otherwise, read dataset into array data[];
 *  - set dataset to corresponding fields in bodies.
 */
template <typename T>
void H5P_bodiesReadDataset(std::vector<body> &bodies, hid_t &file_id,
                           const char *dsname, T *data,
                           size_t dim = IO_nparticlesproc) {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // reset data buffer to zero
  std::fill(data, data + IO_nparticlesproc, 0.);

  // read dataset
  int err = H5P_readDataset(file_id, dsname, data);
  if (err)
    rank || clog(warn) << "Unable to read " << dsname << ": "
                       << "error code " << err << std::endl;

  // assign corresponding field in bodies
  if (!strcmp(dsname, "type")) {
    for (int64_t i = 0; i < IO_nparticlesproc; ++i)
      bodies[i].setType(data[i]);
  } else if (!strcmp(dsname, "id")) {
    if (err == 0) {
      // set existing IDs from file
      for (int64_t i = 0; i < IO_nparticlesproc; ++i)
        bodies[i].set_id(data[i]);
    } else {
      // generate the ids
      rank || clog(trace) << "Setting ID for particles" << std::endl;
      int64_t start = (IO_nparticles / size) * rank + 1;
      for (int64_t i = 0; i < IO_nparticlesproc; ++i) {
        bodies[i].set_id(start + i);
      }
      if (rank == size - 1) // last rank
        assert(IO_nparticles == bodies.back().id());
    }
  } else if (!strcmp(dsname, "m")) {
    for (int64_t i = 0; i < IO_nparticlesproc; ++i)
      bodies[i].set_mass(data[i]);
  } else if (!strcmp(dsname, "rho")) {
    for (int64_t i = 0; i < IO_nparticlesproc; ++i)
      bodies[i].setDensity(data[i]);
  } else if (!strcmp(dsname, "h")) {
    for (int64_t i = 0; i < IO_nparticlesproc; ++i)
      bodies[i].set_radius(data[i]);
  } else if (!strcmp(dsname, "P")) {
    for (int64_t i = 0; i < IO_nparticlesproc; ++i)
      bodies[i].setPressure(data[i]);
  }
#ifdef INTERNAL_ENERGY
  else if (!strcmp(dsname, "u")) {
    for (int64_t i = 0; i < IO_nparticlesproc; ++i)
      bodies[i].setInternalenergy(data[i]);
  }
#endif
  else if (!strcmp(dsname, "dt")) {
    for (int64_t i = 0; i < IO_nparticlesproc; ++i)
      bodies[i].setDt(data[i]);
  } else if constexpr (std::is_same_v<T, double>) {
    if (!strcmp(dsname, "x")) {
      if constexpr (gdimension == 1) {
        for (int64_t i = 0; i < IO_nparticlesproc; ++i) {
          point_t pos = {data[i]};
          bodies[i].set_coordinates(pos);
        }
      }
      if constexpr (gdimension == 2) {
        std::fill(data + IO_nparticlesproc, data + IO_nparticlesproc * 2, 0.);
        H5P_readDataset(file_id, "y", data + IO_nparticlesproc);
        for (int64_t i = 0; i < IO_nparticlesproc; ++i) {
          point_t pos = {data[i], data[IO_nparticlesproc + i]};
          bodies[i].set_coordinates(pos);
        }
      }
      if constexpr (gdimension == 3) {
        std::fill(data + IO_nparticlesproc, data + IO_nparticlesproc * 3, 0.);
        H5P_readDataset(file_id, "y", data + IO_nparticlesproc);
        H5P_readDataset(file_id, "z", data + 2 * IO_nparticlesproc);
        for (int64_t i = 0; i < IO_nparticlesproc; ++i) {
          point_t pos = {data[i], data[IO_nparticlesproc + i],
                         data[IO_nparticlesproc * 2 + i]};
          bodies[i].set_coordinates(pos);
        }
      }
    } else if (!strcmp(dsname, "vx")) {
      if constexpr (gdimension == 1) {
        for (int64_t i = 0; i < IO_nparticlesproc; ++i) {
          point_t vel = {data[i]};
          bodies[i].setVelocity(vel);
        }
      }
      if constexpr (gdimension == 2) {
        std::fill(data + IO_nparticlesproc, data + IO_nparticlesproc * 2, 0.);
        H5P_readDataset(file_id, "vy", data + IO_nparticlesproc);
        for (int64_t i = 0; i < IO_nparticlesproc; ++i) {
          point_t vel = {data[i], data[IO_nparticlesproc + i]};
          bodies[i].setVelocity(vel);
        }
      }
      if constexpr (gdimension == 3) {
        std::fill(data + IO_nparticlesproc, data + IO_nparticlesproc * 3, 0.);
        H5P_readDataset(file_id, "vy", data + IO_nparticlesproc);
        H5P_readDataset(file_id, "vz", data + 2 * IO_nparticlesproc);
        for (int64_t i = 0; i < IO_nparticlesproc; ++i) {
          point_t vel = {data[i], data[IO_nparticlesproc + i],
                         data[IO_nparticlesproc * 2 + i]};
          bodies[i].setVelocity(vel);
        }
      } // switch gdimension
    }   // if dsname
  }     // if T is double

} // H5P_bodiesReadDataset()

size_t H5P_setNumParticles(const int64_t &nparticlesproc) {
  int rank, size;
  MPI_Comm_rank(comm_, &rank);
  MPI_Comm_size(comm_, &size);
  // Compute arrays for hyperslab

  IO_nparticlesproc = nparticlesproc;
  std::vector<int64_t> offcount(size + 1);
  MPI_Allgather(&nparticlesproc, 1, MPI_INT64_T, &offcount[0], 1, MPI_INT64_T,
                comm_);

  for (int i = 1; i < size; ++i)
    offcount[i] += offcount[i - 1];
  offcount.insert(offcount.begin(), 0);

  IO_nparticles = offcount[size];
  IO_offset = offcount[rank];
  IO_count = nparticlesproc;

  return IO_nparticlesproc;
}

size_t H5P_getNumParticles(hid_t file_id) {
  // std::cout<<"Get N PART";
  // Read x and get the dimension of the dataset
  hid_t status = 1;
  // char cstep[255];
  // sprintf(cstep,"/#Step%lu/x",step);
  hid_t dset_id = H5Dopen(IO_group_id, "x", H5P_DEFAULT);
  hid_t dspace = H5Dget_space(dset_id);
  // size_t parts = H5Sget_simple_extent_ndims(dspace);
  hsize_t dim;
  H5Sget_simple_extent_dims(dspace, &dim, NULL);
  size_t parts = dim;
  H5Dclose(dset_id);
  H5Sclose(dspace);
  // std::cout<<"done: "<<parts<<std::endl;
  return parts;
}

/*
 * @brief    Checks whether a filename is a snapshot of the form:
 *           <fprefix>_XXXXX.h5part
 * @param    fprefix    file prefix
 * @param    filename   file name
 * @return   step number
 *
 */
int H5P_isPrefixSnapshot(const char *fprefix, const char *filename) {
  if (strstr(filename, fprefix) == NULL)
    return -1;
  if (strstr(filename, ".h5part") == NULL)
    return -1;
  int imx = strlen(filename) - 8;
  int imn = imx - 5;
  int snum = 0, t = 1;
  for (int i = imx; i > imn; --i) {
    char c = filename[i];
    if (c < '0' or c > '9')
      return -1;
    snum += t * (c - '0');
    t *= 10;
  }
  if (filename[imn] != '_')
    return -1;
  char buf[MAX_FNAME_LEN];
  strcpy(buf, filename);
  buf[imn] = '\0';
  if (strcmp(fprefix, buf) == 0)
    return snum;
  else
    return -1;
}

/*
 * @brief    Checks if file exists, C-style
 * @param    prefix     - filename prefix
 */
bool H5P_fileExists(const char *prefix) {
  char fname[MAX_FNAME_LEN];
  sprintf(fname, "%s.h5part", prefix);
  return (access(fname, F_OK) != -1);
}

/*
 * @brief    For a multiple-files H5part data, returns the snapshot
 *           number of a file with given iteration
 * @param    prefix     - filename prefix
 */
int H5P_findIterationSnapshot(const char *prefix, const int iteration) {
  int step = -1;
  char fname[MAX_FNAME_LEN];
  char prefix_dirname[MAX_FNAME_LEN], buf[MAX_FNAME_LEN],
      prefix_basename[MAX_FNAME_LEN];
  DIR *d;
  struct dirent *dir;
  hid_t file_id = 0;

  // MPI stuff
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // open the prefix directory
  sprintf(buf, prefix);
  sprintf(prefix_dirname, dirname(buf));
  sprintf(buf, prefix);
  sprintf(prefix_basename, basename(buf));
  d = opendir(prefix_dirname);
  if (d) {
    // go through individual files
    while ((dir = readdir(d)) != NULL) {

      // skip files which are not "prefix_XXXXX.h5part"
      step = H5P_isPrefixSnapshot(prefix_basename, dir->d_name);
      if (step < 0)
        continue;

      // if the found file is malformed, abort
      sprintf(fname, "%s_%05d.h5part", prefix, step);
      file_id = H5P_openFile(fname, H5F_ACC_RDONLY);
      if (not H5P_hasStep(file_id, step)) {
        rank || clog(error) << "Cannot find snapshot '" << step << "' in Step#"
                            << step << " in file " << fname << std::endl;
        FULLSTOP;
      }
      H5P_setStep(file_id, step);

      // get iteration of the step
      int64_t file_iteration;
      if (0 != H5P_readAttributeStep(file_id, "iteration", &file_iteration)) {
        rank || clog(error) << "Cannot read attribute 'iteration' in Step#"
                            << step << " in file " << fname << std::endl;
        FULLSTOP;
      }
      H5Fclose(file_id);

      // check if this is what we are looking for
      if (file_iteration == iteration)
        break;
    }
    closedir(d);
  }
  return step;
}

/*
 * @brief    Remove *.h5part files with prefix output_file_prefix
 *           If out_h5data_separate_iterations, remove only steps after
 *           the specified one
 * @param    output_file_prefix    the prefix
 * @param    threshold_stepnum     delete everything > snapshot
 */
int H5P_removePrefix(const char *output_file_prefix,
                     const int threshold_stepnum) {
  char output_filename[MAX_FNAME_LEN];
  int n_deleted = 0;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (not param::out_h5data_separate_iterations) {
    sprintf(output_filename, "%s.h5part", output_file_prefix);
    if (remove(output_filename) == 0) { // if successful
      rank || clog(warn) << "deleting old output file: " << output_filename
                         << std::endl;
      ++n_deleted;
    }
  } else {
    char output_dirname[MAX_FNAME_LEN], buf[MAX_FNAME_LEN],
        output_basename[MAX_FNAME_LEN];
    DIR *d;
    struct dirent *dir;

    sprintf(buf, output_file_prefix);
    sprintf(output_dirname, dirname(buf));
    sprintf(buf, output_file_prefix);
    sprintf(output_basename, basename(buf));
    d = opendir(output_dirname);
    if (d) {
      while ((dir = readdir(d)) != NULL) {
        int stepnum = H5P_isPrefixSnapshot(output_basename, dir->d_name);
        if (stepnum > threshold_stepnum and rank == 0) {
          if (remove(dir->d_name) == 0) { // if successful
            rank || clog(warn) << "deleting old output file: " << dir->d_name
                               << std::endl;
            ++n_deleted;
          }
        }
      }
      closedir(d);
    }
  }
  return n_deleted;
}

// Input data fro HDF5 File
void inputDataHDF5(std::vector<body> &bodies, const char *input_file_prefix,
                   const char *output_file_prefix, int64_t &totalnbodies,
                   int64_t &nbodies, int startIteration,
                   MPI_Comm comm = MPI_COMM_WORLD) {

  comm_ = comm;
  char input_filename[MAX_FNAME_LEN];

  // add the .h5part extension
  sprintf(input_filename, "%s.h5part", input_file_prefix);
  bool input_single_file = H5P_fileExists(input_file_prefix);
  hid_t dataFile;

  // Default if new file, startStep = 0
  int startStep = 0;

  // Lower message level, handle them here on one process
  // H5SetVerbosityLevel  (0);

  int rank, size;
  MPI_Comm_size(comm_, &size);
  MPI_Comm_rank(comm_, &rank);

  clog_one(trace) << "Input particles" << std::endl;

  // ------------- START FROM ITERATION ZERO: OVERWRITE OUTPUT  --------
  if (startIteration == 0) {

    dataFile = H5P_openFile(input_filename, H5F_ACC_RDONLY);

  }
  // ------------- *** OR *** START FROM SPECIFIED ITERATION  ----------
  else {

    // step filename
    char step_filename[MAX_FNAME_LEN];

    if (input_single_file) { // --- single-file mode ---

      // Go through all the steps in a single file
      dataFile = H5P_openFile(input_filename, H5F_ACC_RDONLY);

      // TODO: this loop better run over all timesteps
      const int64_t maxstep = 1000000;
      int64_t step;
      for (step = 0; step < maxstep; ++step) {
        if (H5P_hasStep(dataFile, step)) {
          // Check iteration number
          H5P_setStep(dataFile, step);
          int64_t iteration;
          if (0 != H5P_readAttributeStep(dataFile, "iteration", &iteration)) {
            clog_one(error)
                << "Cannot find attribute 'iteration' in Step#" << step
                << " in file " << input_filename << std::endl;
            FULLSTOP;
          }
          if (iteration == startIteration)
            break;
        }
      }
      if (step == maxstep) {
        clog_one(error) << "Cannot find iteration " << startIteration << " in "
                        << input_filename << std::endl;
        FULLSTOP;
      }
      startStep = step;
      clog_one(info) << "Found iteration " << startIteration << " at step Step#"
                     << startStep << " in " << input_filename << std::endl;
    } else { // ---- multiple-file mode ---

      // find the file with initial_iteration
      int step = H5P_findIterationSnapshot(input_file_prefix,
                                           param::initial_iteration);
      // file doesn't exist: complain and exit
      if (step < 0) {
        clog_one(error) << "Cannot find iteration " << param::initial_iteration
                        << " in prefix " << input_file_prefix << std::endl;
        FULLSTOP;
      }

      // open it
      char step_filename[128];
      sprintf(step_filename, "%s_%05d.h5part", input_file_prefix, step);
      clog_one(warn) << "Reading from file " << step_filename << std::endl;

      // set dataFile and startStep
      dataFile = H5P_openFile(step_filename, H5F_ACC_RDONLY);
      startStep = step;
    }

  } // if specified iteration

  // ------------- CLEAR OUTPUT ------------------------------------------------
  // output prefix is different from input prefix:
  // -> clear output
  if (strcmp(output_file_prefix, input_file_prefix) != 0) {
    H5P_removePrefix(output_file_prefix, -1);
    output_step = 0;
  } else { // --- output prefix == input prefix

    // multiple-files output mode
    if (param::out_h5data_separate_iterations) {

      if (param::initial_iteration == 0) { // invalid input
        clog_one(error)
            << "Invalid combination: cannot have "
            << "initial_iteration = 0 and input prefix == output prefix "
            << "at the same time" << std::endl;
        FULLSTOP
      } else {
        // only remove files with output prefix that have step numbers higher
        // than then one with intial_iteration
        H5P_removePrefix(output_file_prefix, startStep);
      }
    }
    // single-file output mode
    else {

      // check if the lastStep is the startIteration
      char output_filename[MAX_FNAME_LEN];
      sprintf(output_filename, "%s.h5part", output_file_prefix);
      hid_t outputFile = H5P_openFile(output_filename, H5F_ACC_RDONLY);
      H5P_setStep(outputFile, startStep);
      // Check if startStep == lastStep
      int lastStep = startStep;
      while (H5P_hasStep(outputFile, lastStep)) {
        lastStep++;
      }
      --lastStep;
      clog_one(trace) << "startStep: " << startStep << " lastStep: " << lastStep
                      << std::endl;
      if (startStep != lastStep) {
        clog_one(error) << "First step not last step in output" << std::endl;
        H5P_closeFile(outputFile);
        MPI_Barrier(comm_);
        MPI_Finalize();
      }
      H5P_closeFile(outputFile);
    }
    output_step = startStep;
  } // if input_prefix == output_prefix

  if (dataFile == 0) {
    clog_one(error) << "Cannot find data file" << std::endl;
    MPI_Barrier(comm_);
    MPI_Finalize();
  }

  int val = H5P_hasStep(dataFile, startStep);

  // If yes, get into this step
  H5P_setStep(dataFile, startStep);
  // Get the number of particles
  int64_t nparticles = 0UL;
  nparticles = H5P_getNumParticles(dataFile);

  int64_t nparticlesproc = nparticles / size;
  // Handle the number of particles for the last one
  if (size == rank + 1) {
    nparticlesproc = nparticles - nparticlesproc * (size - 1);
  }

  H5P_setNumParticles(nparticlesproc);

  // Register for main
  totalnbodies = nparticles;
  nbodies = IO_nparticlesproc;

  //--------------- READ GLOBAL ATTRIBUTES ------------------------------------
  // read the number of dimension
  int32_t dimension;
  if (0 == H5P_readAttribute(dataFile, "dimension", &dimension))
    assert(gdimension == dimension);

  //--------------- READ DATA FROM STEP ---------------------------------------
  // Set the number of particles read by each process
  // H5PartSetNumParticles(dataFile,nparticlesproc);
  // Resize the bodies vector
  bodies.clear();
  bodies.resize(IO_nparticlesproc);

  // Try to read timestep if exists
  if (startIteration != 0) {
    double timestep = 1.;
    double totaltime = 0.;
    if (0 == H5P_readAttributeStep(dataFile, "timestep", &timestep)) {
      physics::dt = timestep;
    } else {
      clog_one(warn) << "Attribute 'timestep' missing in input file"
                     << input_filename << std::endl;
    }
    if (0 == H5P_readAttributeStep(dataFile, "time", &totaltime)) {
      physics::totaltime = totaltime;
    } else {
      clog_one(warn) << "Attribute 'totaltime' missing in input file"
                     << input_filename << std::endl;
    }
  }

  // Read the dataset and fill the particles data
  double *dataX = new double[IO_nparticlesproc * gdimension];
  int64_t *dataInt = new int64_t[IO_nparticlesproc];
  int *dataInt32 = new int[IO_nparticlesproc];

  // Read positions and velocities
  H5P_bodiesReadDataset(bodies, dataFile, "x", dataX);
  H5P_bodiesReadDataset(bodies, dataFile, "vx", dataX);
  H5P_bodiesReadDataset(bodies, dataFile, "m", dataX);
  H5P_bodiesReadDataset(bodies, dataFile, "rho", dataX);
  H5P_bodiesReadDataset(bodies, dataFile, "h", dataX);
  H5P_bodiesReadDataset(bodies, dataFile, "P", dataX);

#ifdef INTERNAL_ENERGY
  H5P_bodiesReadDataset(bodies, dataFile, "u", dataX);
#endif

  H5P_bodiesReadDataset(bodies, dataFile, "id", dataInt);
  H5P_bodiesReadDataset(bodies, dataFile, "dt", dataX);
  H5P_bodiesReadDataset(bodies, dataFile, "type", dataInt32);

  delete[] dataX;
  delete[] dataInt;
  delete[] dataInt32;

  MPI_Barrier(comm_);
  H5Fclose(dataFile);
  // H5CloseFile(dataFile);

  clog_one(trace) << "Input particles.done" << std::endl;

} // inputDataHDF5

// Output data in HDF5 format
// Generate the associate XDMF file
void outputDataHDF5(std::vector<body> &bodies, const char *fileprefix,
                    int64_t iteration, double totaltime,
                    MPI_Comm comm = MPI_COMM_WORLD) {

  comm_ = comm;
  int step = output_step++;

  int size, rank;
  MPI_Comm_size(comm_, &size);
  MPI_Comm_rank(comm_, &rank);

  MPI_Barrier(comm_);
  clog_one(trace) << "Output particles" << std::flush;

  char filename[128];
  if (param::out_h5data_separate_iterations)
    sprintf(filename, "%s_%05d.h5part", fileprefix, step);
  else
    sprintf(filename, "%s.h5part", fileprefix);

  // Wait for removing the file before writing in
  MPI_Barrier(comm_);
  // Check if file exists
  hid_t dataFile = H5P_openFile(filename, H5F_ACC_RDWR);

  //-------------------GLOBAL HEADER-------------------------------------------
  // Only for the first output
  if (step == 0 or param::out_h5data_separate_iterations) {
    int gdimension32 = gdimension;
    H5P_writeAttribute(dataFile, "dimension", &gdimension32);
  }

  //------------------STEP HEADER----------------------------------------------
  // Put the step header
  H5P_setStep(dataFile, step);
  H5P_writeAttributeStep(dataFile, "time", &physics::totaltime);
  H5P_writeAttributeStep(dataFile, "iteration", &iteration);
  H5P_writeAttributeStep(dataFile, "timestep", &physics::dt);
  //------------------STEP DATA------------------------------------------------

  H5P_setNumParticles(bodies.size());

  // 4 buffers, 3 double and 1 int64
  double *b1 = new double[IO_nparticlesproc];
  double *b2 = new double[IO_nparticlesproc];
  double *b3 = new double[IO_nparticlesproc];
  int64_t *bi = new int64_t[IO_nparticlesproc];
  int32_t *bint = new int32_t[IO_nparticlesproc];

  // Position
  int64_t pos = 0L;
  // Extract data from bodies
  for (auto bi : bodies) {
    b1[pos] = bi.coordinates()[0];
    if (gdimension > 1) {
      b2[pos] = bi.coordinates()[1];
    } else {
      b2[pos] = 0.;
    }
    if (gdimension > 2) {
      b3[pos++] = bi.coordinates()[2];
    } else {
      b3[pos++] = 0.;
    }
  }
  H5P_writeDataset(dataFile, "x", b1);
  H5P_writeDataset(dataFile, "y", b2);
  H5P_writeDataset(dataFile, "z", b3);

  // Velocity
  pos = 0L;
  // Extract data from bodies
  for (auto bi : bodies) {
    b1[pos] = bi.getVelocity()[0];
    if (gdimension > 1) {
      b2[pos] = bi.getVelocity()[1];
    } else {
      b2[pos] = 0;
    }
    if (gdimension > 2) {
      b3[pos++] = bi.getVelocity()[2];
    } else {
      b3[pos++] = 0.;
    }
  }
  H5P_writeDataset(dataFile, "vx", b1);
  H5P_writeDataset(dataFile, "vy", b2);
  H5P_writeDataset(dataFile, "vz", b3);

  // Acceleration
  pos = 0L;
  // Extract data from bodies
  for (auto bi : bodies) {
    b1[pos] = bi.getAcceleration()[0];
    if (gdimension > 1) {
      b2[pos] = bi.getAcceleration()[1];
    } else {
      b2[pos] = 0.;
    }
    if (gdimension > 2) {
      b3[pos++] = bi.getAcceleration()[2];
    } else {
      b3[pos++] = 0.;
    }
  }
  H5P_writeDataset(dataFile, "ax", b1);
  H5P_writeDataset(dataFile, "ay", b2);
  H5P_writeDataset(dataFile, "az", b3);

  // Smoothing length, Density, Internal Energy
  pos = 0L;
  // Extract data from bodies
  for (auto bi : bodies) {
    b1[pos] = bi.radius();
    b2[pos] = bi.getDensity();
#ifdef INTERNAL_ENERGY
    b3[pos] = bi.getInternalenergy();
#endif
    pos++;
  }
  H5P_writeDataset(dataFile, "h", b1);
  H5P_writeDataset(dataFile, "rho", b2);
#ifdef INTERNAL_ENERGY
  H5P_writeDataset(dataFile, "u", b3);
#endif

  // Pressure, Mass, Id, timestep
  pos = 0L;
  // Extract data from bodies
  for (auto bid : bodies) {
    b1[pos] = bid.getPressure();
    b2[pos] = bid.mass();
    b3[pos] = bid.getDt();
    bi[pos] = bid.id();
    bint[pos++] = bid.getType();
  }
  H5P_writeDataset(dataFile, "P", b1);
  H5P_writeDataset(dataFile, "m", b2);
  H5P_writeDataset(dataFile, "dt", b3);
  H5P_writeDataset(dataFile, "id", bi);
  H5P_writeDataset(dataFile, "type", bint);

  // Output the rank for analysis
  std::fill(bi, bi + IO_nparticlesproc, rank);
  H5P_writeDataset(dataFile, "rank", bi);

  pos = 0L;
  for (auto bid : bodies) {
    bi[pos++] = bid.key().value();
  }
  H5P_writeDataset(dataFile, "key", bi);

  H5P_closeFile(dataFile);

  delete[] b1;
  delete[] b2;
  delete[] b3;
  delete[] bi;
  delete[] bint;

  clog_one(trace) << ".done" << std::endl;

} // outputDataHDF5

} // namespace io

#undef FULLSTOP
#endif // _mpisph_io_h_
