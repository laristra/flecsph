/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>

#include "user.h"
#include "sodtube.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"
#include "io.h"
using namespace io;

//
// help message
//
void print_usage() {
  clog_one(warn)
      << "Initial data generator for Sod shocktube test in"
      << gdimension << "D" << std::endl
      << "Usage: ./sodtube_generator <parameter-file.par>" << std::endl;
}

//
// derived parameters
//
static int64_t nparticlesproc;        // number of particles per proc
static double rho_1, rho_2;           // densities
static double vx_1, vx_2;             // velocities
static double pressure_1, pressure_2; // pressures
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

// geometric extents of the three regions: left, right and central
static point_t cbox_min, cbox_max;
static point_t rbox_min, rbox_max;
static point_t lbox_min, lbox_max;

void set_derived_params() {
  using namespace std;
  using namespace param;

  // test selector
  switch (sodtest_num) {
    case (1):
      // -- middle         | left and right side -- //
      rho_1      = 1.0;      rho_2      = 0.125;
      pressure_1 = 1.0;      pressure_2 = 0.1;
      vx_1       = 0.0;      vx_2       = 0.0;
      break;

    case (2):
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 0.4;      pressure_2 = 0.4;
      vx_1       =-2.0;      vx_2       = 2.0;
      break;

    case (3):
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 1000.;    pressure_2 = 0.01;
      vx_1       = 0.0;      vx_2       = 0.0;
      break;

    case (4):
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 0.01;     pressure_2 = 100.;
      vx_1       = 0.0;      vx_2       = 0.0;
      break;

    case (5):
      rho_1      = 5.99924;  rho_2      = 5.99242;
      pressure_1 = 460.894;  pressure_2 = 46.0950;
      vx_1       = 19.5975;  vx_2       =-6.19633;
      break;

    case (6): // 1D equivalent to the Noh problem
      rho_1      = 1.0;      rho_2      = 1.0;
      pressure_1 = 1.e-6;    pressure_2 = 1.e-6;
      vx_1       = 1.0;      vx_2       =-1.0;
      break;

    default:
      clog_one(error) << "ERROR: invalid test (" << sodtest_num << ")." << endl;
      MPI_Finalize();
      exit(-1);

  }

  // adjust lattice_nx such that it gives 1 in remainder if divided by 3
  SET_PARAM(lattice_nx, ((lattice_nx - 1)/3) * 3 + 1);
  
  // particle spacing
  SET_PARAM(sph_separation, (box_length/(double)(lattice_nx - 1)));

  // compute the total number of particles
  int64_t npd = lattice_nx;

  // 1D setup
  cbox_max[0] =  box_length/6.;
  cbox_min[0] = -box_length/6.;
  lbox_min[0] = -box_length/2.;
  lbox_max[0] =  cbox_min[0];
  rbox_min[0] =  cbox_max[0];
  rbox_max[0] = -lbox_min[0];

  // 2D case
  if (gdimension>1) {
     cbox_max[1] = lbox_max[1] = rbox_max[1] = box_width/2.0;
     cbox_min[1] = lbox_min[1] = rbox_min[1] =-box_width/2.0;
     npd *= (int64_t)((double)lattice_nx*box_width/box_length);
  }

  // 3D case
  if (gdimension>2) {
     cbox_max[2] = lbox_max[2] = rbox_max[2] = box_height/2.0;
     cbox_min[2] = lbox_min[2] = rbox_min[2] =-box_height/2.0;
     npd *= (int64_t)((double)lattice_nx*box_height/box_length);
  }
  SET_PARAM(nparticles, npd);

  // warn about not using periodic boundaries
  if constexpr (gdimension == 2) {
    if (not (periodic_boundary_x and periodic_boundary_y))
      clog_one(warn)
        << "This test is best done with periodic boundaries. Make sure to "
        << " set periodic_boundary_x = yes and periodic_boundary_y = yes"
        << " for this test" << std::endl;
  }
  if constexpr (gdimension == 3) {
    if (not (periodic_boundary_x and periodic_boundary_y and periodic_boundary_z))
      clog_one(warn)
        << "This test is best done with periodic boundaries. Make sure to set:"
        << std::endl << "  periodic_boundary_x = yes"
        << std::endl << "  periodic_boundary_y = yes"
        << std::endl << "  periodic_boundary_z = yes"
        << std::endl << "for best results." << std::endl;
  }

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace param;
  const double b_tol = particle_lattice::b_tol;

  // launch MPI
  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // check options list: exactly one option is allowed
  if (argc != 2) {
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // screen output
  clog_one(info)
    << "Sod shocktube test #" << sodtest_num
    << "in " << gdimension << "D" << std::endl;

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();
  particle_lattice::select();

  // set kernel
  kernels::select();

  // allocate arrays
  int64_t tparticles = 0;
  int64_t parts_mid = 0;
  int64_t parts_lr = 0;
  double mass = 0;
  bool equal_separation = !equal_mass;
  if(equal_separation){
    tparticles =  particle_lattice::count(lattice_type,2,cbox_min,cbox_max,
                                          sph_separation,0);
    parts_mid = tparticles;
    tparticles += particle_lattice::count(lattice_type,2,rbox_min,rbox_max,
                                          sph_separation,tparticles-1);
    parts_lr = tparticles-parts_mid;
    tparticles += particle_lattice::count(lattice_type,2,lbox_min,lbox_max,
                                          sph_separation,tparticles-1);
  }

  // screen output
  clog_one(info)
    << "Number of particles: "
    << tparticles << std::endl;
  clog_one(info)
    << "Initial data file: " 
    << initial_data_file << std::endl;

  double lr_sph_sep = 0.;
  double temp_part = 0;
  double temp_part_new = 0;
  double dy1, dy2, dz1, dz2;
  dy1 = dy2 = dz1 = dz2 = sph_separation;
  if(equal_mass){
    if (gdimension == 1) {
      mass = rho_1*sph_separation;
      lr_sph_sep = mass/rho_2;
      clog_one(info) << std::endl 
        << "Lattice resolution: " << std::endl
        << " - central box:      dx = " << sph_separation << std::endl
        << " - left/right boxes: dx = " << lr_sph_sep << std::endl;
    } 
    else if (gdimension == 2) {
      mass = rho_1*sph_separation*sph_separation;
      if (lattice_type == 1 or lattice_type == 2)
        mass *= sqrt(3.0)/2.0;
      lr_sph_sep = sph_separation * sqrt(rho_1/rho_2);
      if (lattice_type == 0) {
        dy1 = sph_separation;
        dy2 = lr_sph_sep;
        clog_one(info) << std::endl 
          << "Lattice resolution: " << std::endl
          << " - central box:      dx = " << sph_separation << std::endl
          << " - left/right boxes: dx = " << lr_sph_sep << std::endl;
      }
      else if (lattice_type == 1 or lattice_type == 2) {
        dy1 = sph_separation*sqrt(3.);
        dy2 = lr_sph_sep*sqrt(3.);
        clog_one(info) << std::endl 
          << "Lattice resolution: " << std::endl
          << " - central box:      dx = " << sph_separation << std::endl
          << "                   2*dy = " << dy1 << std::endl
          << " - left/right boxes: dx = " << lr_sph_sep << std::endl
          << "                   2*dy = " << dy2 << std::endl;
      }

    } 
    else { // dimension == 3
      mass = rho_1*sph_separation*sph_separation*sph_separation;
      if (lattice_type == 1 or lattice_type == 2)
        mass *= 1.0/sqrt(2.0);
      lr_sph_sep = sph_separation * cbrt(rho_1/rho_2);
      if (lattice_type == 0) {
        dy1 = dz1 = sph_separation;
        dy2 = dz2 = lr_sph_sep;
        clog_one(info) << std::endl 
          << "Lattice: rectangular, resolution: " << std::endl
          << " - central box:      dx = " << sph_separation << std::endl
          << " - left/right boxes: dx = " << lr_sph_sep << std::endl;
      }
      else if (lattice_type == 1) {
        dy1 = sph_separation*sqrt(3.);
        dy2 = lr_sph_sep*sqrt(3.);
        dz1 = 2.*sph_separation*sqrt(2./3.);
        dz2 = 2.*lr_sph_sep*sqrt(2./3.);
        clog_one(info) << std::endl 
          << "Lattice: HCP, resolution: " << std::endl
          << " - central box:      dx = " << sph_separation << std::endl
          << "                   2*dy = " << dy1 << std::endl
          << "                   2*dz = " << dz1 << std::endl
          << " - left/right boxes: dx = " << lr_sph_sep << std::endl
          << "                   2*dy = " << dy2 << std::endl
          << "                   2*dz = " << dz2 << std::endl;
      }
      else if (lattice_type == 2) {
        dy1 = sph_separation*sqrt(3.);
        dy2 = lr_sph_sep*sqrt(3.);
        dz1 = sph_separation*sqrt(6.);
        dz2 = lr_sph_sep*sqrt(6.);
        clog_one(info) << std::endl 
          << "Lattice: FCC, resolution: " << std::endl
          << " - central box:      dx = " << sph_separation << std::endl
          << "                   2*dy = " << dy1 << std::endl
          << "                   3*dz = " << dz1 << std::endl
          << " - left/right boxes: dx = " << lr_sph_sep << std::endl
          << "                   2*dy = " << dy2 << std::endl
          << "                   3*dz = " << dz2 << std::endl;
      }
    }

    // for periodic boundaries, lattice has to match up:
    // adjust domain length
    if (periodic_boundary_x) {
      SET_PARAM(box_length, 
                box_length/3 + 2*(int(box_length/(3*lr_sph_sep)))*lr_sph_sep)
      rbox_max[0] = box_length/2;
      lbox_min[0] =-box_length/2;
    }
    
    // adjust domain width
    if (gdimension >= 2 and periodic_boundary_y) {
      int Ny1 = (int)(box_width/dy1) - 1;
      for(int j = Ny1; j < Ny1*100; ++j) {
        double w2 = (floor(j*dy1/dy2 - b_tol) + 1)*dy2;
        if (fabs(w2 - j*dy1) < lattice_mismatch_tolerance*dy2) {
          SET_PARAM(box_width, std::min(w2,j*dy1));
          cbox_min[1] = -box_width/2.;
          cbox_max[1] =  box_width/2.;
          rbox_min[1] = -box_width/2.;
          rbox_max[1] =  box_width/2.;
          lbox_min[1] = -box_width/2.;
          lbox_max[1] =  box_width/2.;
          break;
        }
      }
    }

    // adjust domain height
    if (gdimension >= 3 and periodic_boundary_z) {
      int Nz1 = (int)(box_height/dz1) - 1;
      for(int k = Nz1; k < Nz1*100; ++k) {
        double w2 = (floor(k*dz1/dz2 - b_tol) + 1)*dz2;
        if (fabs(w2 - k*dz1) < lattice_mismatch_tolerance*dz2) {
          SET_PARAM(box_height, std::min(w2,k*dz1));
          cbox_min[2] = -box_height/2.;
          cbox_max[2] =  box_height/2.;
          rbox_min[2] = -box_height/2.;
          rbox_max[2] =  box_height/2.;
          lbox_min[2] = -box_height/2.;
          lbox_max[2] =  box_height/2.;
          break;
        }
      }
    }

    // report adjusted dimensions
    if (periodic_boundary_x or periodic_boundary_y or periodic_boundary_z) {
      clog_one(warn) 
        << "Domain has been adjusted for periodic boundaries." << std::endl;
      clog_one(warn) 
        << "For evolution, modify domain dimensions as follows:" << std::endl 
        << "  box_length = " << box_length << std::endl 
        << "  box_width = "  << box_width  << std::endl 
        << "  box_height = " << box_height << std::endl;
      
      //if constexpr (gdimension >= 1) 
      //  clog_one(warn) << "Lattice mismatch on the boundaries:" << std::endl
      //    << " X-direction: " 
      //    << (box_length - cbox_min[0] + lbox_max[0] 
      //                   - rbox_min[0] + cbox_max[0]
      //     - floor((cbox_max[0] - cbox_min[0])/sph_separation)*sph_separation
      //     - floor((lbox_max[0] - lbox_min[0])/lr_sph_sep)*lr_sph_sep
      //     - floor((rbox_max[0] - rbox_min[0])/lr_sph_sep)*lr_sph_sep)
      //    <<   ", dx = " << sph_separation 
      //    << std::endl;

      if constexpr (gdimension >= 2) 
        clog_one(warn) << "Lattice mismatch, Y-direction:" << std::endl
          << " - central box: "
          <<   (box_width-floor((cbox_max[1]-cbox_min[1])/dy1)*dy1) 
          <<   ", dy = " << dy1 << ", mismatch/dy = " 
          <<   (box_width/dy1-floor((cbox_max[1]-cbox_min[1])/dy1)) 
          << std::endl
          << " -    left box: "
          <<   (box_width-floor((lbox_max[1]-lbox_min[1])/dy2)*dy2) 
          <<   ", dy = " << dy2 << ", mismatch/dy = " 
          <<   (box_width/dy2-floor((lbox_max[1]-lbox_min[1])/dy2)) 
          << std::endl
          << " -   right box: "
          <<   (box_width-floor((rbox_max[1]-rbox_min[1])/dy2)*dy2)
          <<   ", dy = " << dy2 << ", mismatch/dy = " 
          <<   (box_width/dy2-floor((rbox_max[1]-rbox_min[1])/dy2)) 
          << std::endl;

      if constexpr (gdimension >= 3) 
        clog_one(warn) << "Lattice mismatch, Z-direction:" << std::endl
          << " - central box: "
          <<   (box_height-floor((cbox_max[2]-cbox_min[2])/dz1)*dz1) 
          <<   ", dz = " << dz1 << ", mismatch/dz = " 
          <<   (box_height/dz1-floor((cbox_max[2]-cbox_min[2])/dz1)) 
          << std::endl
          << " -    left box: "
          <<   (box_height-floor((lbox_max[2]-lbox_min[2])/dz2)*dz2) 
          <<   ", dz = " << dz2 << ", mismatch/dz = " 
          <<   (box_height/dz2-floor((lbox_max[2]-lbox_min[2])/dz2)) 
          << std::endl
          << " -   right box: "
          <<   (box_height-floor((rbox_max[2]-rbox_min[2])/dz2)*dz2) 
          <<   ", dz = " << dz2 << ", mismatch/dz = " 
          <<   (box_height/dz2-floor((rbox_max[2]-rbox_min[2])/dz2)) 
          << std::endl;

    }

    tparticles = particle_lattice::count(lattice_type,2,cbox_min,cbox_max,
                                         sph_separation,0);
    tparticles += particle_lattice::count(lattice_type,2,rbox_min,rbox_max,
                                          lr_sph_sep,tparticles);
    tparticles += particle_lattice::count(lattice_type,2,lbox_min,lbox_max,
                                          lr_sph_sep,tparticles);
  } // equal mass

  // Initialize the arrays to be filled later
  // Position
  double* x = new double[tparticles]();
  double* y = new double[tparticles]();
  double* z = new double[tparticles]();
  // Velocity
  double* vx = new double[tparticles]();
  double* vy = new double[tparticles]();
  double* vz = new double[tparticles]();
  // Acceleration
  double* ax = new double[tparticles]();
  double* ay = new double[tparticles]();
  double* az = new double[tparticles]();
  // Smoothing length
  double* h = new double[tparticles]();
  // Density
  double* rho = new double[tparticles]();
  // Internal Energy
  double* u = new double[tparticles]();
  // Pressure
  double* P = new double[tparticles]();
  // Mass
  double* m = new double[tparticles]();
  // Id
  int64_t* id = new int64_t[tparticles]();
  // Timestep
  double* dt = new double[tparticles]();

  if(equal_separation){
    tparticles =  particle_lattice::generate(lattice_type,2,cbox_min,cbox_max,
                                             sph_separation,0,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,rbox_min,rbox_max,
                                             sph_separation,tparticles,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,lbox_min,lbox_max,
                                             sph_separation,tparticles,x,y,z);
  }
  if(equal_mass){
    tparticles =  particle_lattice::generate(lattice_type,2,cbox_min,cbox_max,
                                             sph_separation,0,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,rbox_min,rbox_max,
                                             lr_sph_sep,tparticles,x,y,z);
    tparticles += particle_lattice::generate(lattice_type,2,lbox_min,lbox_max,
                                             lr_sph_sep,tparticles,x,y,z);
  }

  // particle id number
  int64_t posid = 0;

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*std::max(pressure_1/rho_1,pressure_2/rho_2));

  // The value for constant timestep
  double timestep = 0.5*sph_separation/cs;

  if(equal_separation){
    for(int64_t part=0; part<tparticles; ++part){
      id[part] = posid++;
      if (particle_lattice::in_domain_1d(x[part],
          cbox_min[0], cbox_max[0], domain_type)) {
        P[part] = pressure_1;
        rho[part] = rho_1;
        vx[part] = vx_1;
        m[part] = rho[part]/(double)parts_mid;
      }
      else {
        P[part] = pressure_2;
        rho[part] = rho_2;
        vx[part] = vx_2;
        m[part] = rho[part]/(double)parts_lr;
      }

      // compute internal energy using gamma-law eos
      u[part] = P[part]/(poly_gamma-1.)/rho[part];

      // particle smoothing length
      h[part] = sph_eta * kernels::kernel_width
                        * pow(m[part]/rho[part],1./gdimension);
    } // for part=0..nparticles
  }
  if(equal_mass){
    for(int64_t part=0; part<tparticles; ++part){
      id[part] = posid++;
      if (particle_lattice::in_domain_1d(x[part],
          cbox_min[0], cbox_max[0], domain_type)) {
        P[part] = pressure_1;
        rho[part] = rho_1;
        vx[part] = vx_1;
      }
      else {
        P[part] = pressure_2;
        rho[part] = rho_2;
        vx[part] = vx_2;
      }

      // compute internal energy using gamma-law eos
      u[part] = P[part]/(poly_gamma-1.)/rho[part];

      // particle masses and smoothing length
      m[part] = mass;
      h[part] = sph_eta * kernels::kernel_width
                        * pow(m[part]/rho[part],1./gdimension);

    } // for part=0..nparticles
  }
  clog_one(info) << "Actual number of particles: " << tparticles << std::endl
    << std::flush;
  // delete the output file if exists
  remove(initial_data_file.c_str());
  hid_t dataFile = H5P_openFile(initial_data_file.c_str(),H5F_ACC_RDWR);

  int use_fixed_timestep = 1;
  // add the global attributes
  H5P_writeAttribute(dataFile,"nparticles",&tparticles);
  H5P_writeAttribute(dataFile,"timestep",&timestep);
  int dim = gdimension;
  H5P_writeAttribute(dataFile,"dimension",&dim);
  H5P_writeAttribute(dataFile,"use_fixed_timestep",&use_fixed_timestep);

  H5P_setNumParticles(tparticles);
  H5P_setStep(dataFile,0);

  //H5PartSetNumParticles(dataFile,nparticles);
  H5P_writeDataset(dataFile,"x",x,tparticles);
  H5P_writeDataset(dataFile,"y",y,tparticles);
  H5P_writeDataset(dataFile,"z",z,tparticles);
  H5P_writeDataset(dataFile,"vx",vx,tparticles);
  H5P_writeDataset(dataFile,"vy",vy,tparticles);
  H5P_writeDataset(dataFile,"h",h,tparticles);
  H5P_writeDataset(dataFile,"rho",rho,tparticles);
  H5P_writeDataset(dataFile,"u",u,tparticles);
  H5P_writeDataset(dataFile,"P",P,tparticles);
  H5P_writeDataset(dataFile,"m",m,tparticles);
  H5P_writeDataset(dataFile,"id",id,tparticles);

  H5P_closeFile(dataFile);

  delete[] x, y, z, vx, vy, vz, ax, ay, az, h, rho, u, P, m, id, dt;

  MPI_Finalize();
  return 0;
}
