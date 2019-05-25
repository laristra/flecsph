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
      << "Initial data generator for KH test in"
      << gdimension << "D" << std::endl
      << "Usage: ./KD_XD_generator <parameter-file.par>" << std::endl;
}

//
// derived parameters
//
static int64_t nparticlesproc;        // number of particles per proc
static double rho_m, rho_t;           // densities
static double vx_m, vx_t;             // velocities
static double pressure_m, pressure_t; // pressures
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

// geometric extents of the three regions: top, middle and bottom
static point_t tbox_min, tbox_max;
static point_t mbox_min, mbox_max;
static point_t bbox_min, bbox_max;

static int64_t np_middle = 0; // number of particles in the middle block
static int64_t np_top = 0;    // number of particles in the top block
static int64_t np_bottom = 0; // number of particles in the bottom block
static double sph_sep_t = 0; // particle separation in top or bottom blocks
static double pmass = 0;      // particle mass in the middle block
static double pmass_t = 0;   // particle mass in top or bottom blocks

void set_derived_params() {
  using namespace std;
  using namespace param;

  // boundary tolerance factor
  const double b_tol = particle_lattice::b_tol;

  // support for only equal-mass configurations for now
  if(not equal_mass) {
    clog_one(error) 
      << "Only equal-mass configurations are implemented"
      << std::endl;                    
    MPI_Finalize();
    exit(0);
  }

  mbox_max[0] = bbox_max[0] = tbox_max[0] = box_length/2.;
  mbox_min[0] = bbox_min[0] = tbox_min[0] =-box_length/2.;

  mbox_max[1] =  box_width/4.;
  mbox_min[1] = -box_width/4.;
  bbox_min[1] = -box_width/2.;
  bbox_max[1] = -box_width/4.;
  tbox_min[1] =  box_width/4.;
  tbox_max[1] =  box_width/2.;

  if(gdimension == 3){
    mbox_max[2] = bbox_max[2] = tbox_max[2] = box_height/2.;
    mbox_min[2] = bbox_min[2] = tbox_min[2] =-box_height/2.;
  }


  // set physical parameters
  // --  in the top and bottom boxes (tbox and bbox):
  rho_t      =  rho_initial;
  pressure_t =  pressure_initial;
  vx_t       = -flow_velocity/2.0;
  // -- in the middle box
  rho_m      =  rho_t * density_ratio; // 2.0 by default
  pressure_m =  pressure_t;          // pressures must be equal in KH test
  vx_m       =  flow_velocity/2.0;

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

  // select particle lattice and kernel function
  particle_lattice::select();
  kernels::select();

  // particle mass and spacing
  SET_PARAM(sph_separation, box_length/(double)(lattice_nx-1));
  if(gdimension == 3){
    pmass = rho_m*sph_separation*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= 1./sqrt(2.);
    sph_sep_t = sph_separation*cbrt(rho_m/rho_t);
  }
  if(gdimension == 2){
    pmass = rho_m*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= sqrt(3)/2;
    sph_sep_t = sph_separation*sqrt(rho_m/rho_t);
  }
  pmass_t = pmass;

  // lattice spacing
  double dx_m, dy_m, dz_m, dx_t, dy_t, dz_t;
  dx_m = dy_m = dz_m = sph_separation;
  dx_t = dy_t = dz_t = sph_sep_t;
  if (lattice_type == 0) {
    clog_one(info)
      << "Lattice: rectangular, resolution: " << std::endl
      << " - middle box:     dx = " << dx_m << std::endl
      << " - top/bottom box: dx = " << dx_t << std::endl;
  }
  if (lattice_type == 1) { // HCP lattice
    dy_m *= sqrt(3.);
    dy_t *= sqrt(3.);
    dz_m *= 2.*sqrt(2./3.);
    dz_t *= 2.*sqrt(2./3.);
    clog_one(info)
      << "Lattice: HCP, resolution: " << std::endl
      << " - middle box:     dx = " << dx_m << std::endl
      << "                 2*dy = " << dy_m << std::endl
      << "                 2*dz = " << dz_m << std::endl
      << " - top/bottom box: dx = " << dx_t << std::endl
      << "                 2*dy = " << dy_t << std::endl
      << "                 2*dz = " << dz_t << std::endl;
  }
  if (lattice_type == 2) { // FCC lattice
    dy_m *= sqrt(3.);
    dy_t *= sqrt(3.);
    dz_m *= 3.*sqrt(2./3.);
    dz_t *= 3.*sqrt(2./3.);
    clog_one(info)
      << "Lattice: FCC, resolution: " << std::endl
      << " - middle box:     dx = " << dx_m << std::endl
      << "                 2*dy = " << dy_m << std::endl
      << "                 3*dz = " << dz_m << std::endl
      << " - top/bottom box: dx = " << dx_t << std::endl
      << "                 2*dy = " << dy_t << std::endl
      << "                 3*dz = " << dz_t << std::endl;
  }

  // adjust width in y-direction of the middle block for symmetry
  double w_m = ((int)(box_width/(3.*dy_m)))*dy_m;
  mbox_min[1] = -0.5*w_m;
  mbox_max[1] =  0.5*w_m + 2*b_tol;

  // adjust top and bottom blocks
  double w_t = ((int)((box_width/2. - w_m/2.)/dy_t))*dy_t;
  double gap = std::min(dy_m,dy_t)/2.;
  SET_PARAM(box_width, w_m + 2.*(w_t + gap))
  tbox_min[1] = 0.5*box_width - w_t;
  tbox_max[1] = 0.5*box_width;
  bbox_min[1] =-0.5*box_width;
  bbox_max[1] =-0.5*box_width + w_t + 2*b_tol;

  // adjust domain length
  int Nx = (int)(box_length/dx_m) - 1;
  for(int i = Nx; i < Nx*100; ++i) {
    double w2 = floor(i*dx_m/dx_t + b_tol)*dx_t;
    if (fabs(w2 - i*dx_m) < lattice_mismatch_tolerance*dx_t) {
      SET_PARAM(box_length, std::min(w2,i*dx_m));
      bbox_min[0] = -box_length/2.;
      bbox_max[0] =  box_length/2.;
      mbox_min[0] = -box_length/2.;
      mbox_max[0] =  box_length/2.;
      tbox_min[0] = -box_length/2.;
      tbox_max[0] =  box_length/2.;
      break;
    }
  }

  // adjust domain height
  if constexpr (gdimension >= 3) {
    int Nz = (int)(box_height/dz_m) - 1;
    for(int k = Nz; k < Nz*100; ++k) {
      double w2 = floor(k*dz_m/dz_t + b_tol)*dz_t;
      if (fabs(w2 - k*dz_m) < lattice_mismatch_tolerance*dz_t) {
        SET_PARAM(box_height, std::min(w2,k*dz_m));
        bbox_min[2] = -box_height/2.;
        bbox_max[2] =  box_height/2.;
        mbox_min[2] = -box_height/2.;
        mbox_max[2] =  box_height/2.;
        tbox_min[2] = -box_height/2.;
        tbox_max[2] =  box_height/2.;
        break;
      }
    }
  }

  // report adjusted dimensions
  clog_one(warn) 
    << "Domain has been adjusted for periodic boundaries." << std::endl;
  clog_one(warn) 
    << "For evolution, modify domain dimensions as follows:" << std::endl 
    << "  box_length = " << box_length << std::endl 
    << "  box_width = "  << box_width  << std::endl 
    << "  box_height = " << box_height << std::endl;


  clog_one(warn) << "Lattice mismatch, X-direction:" << std::endl
      << " -    top box: "
      <<   (box_length - floor((tbox_max[0]-tbox_min[0])/dx_t)*dx_t) 
      <<   ", dx = " << dx_t << ", mismatch/dx = " 
      <<   (box_length/dx_t - floor((tbox_max[0]-tbox_min[0])/dx_t)) 
      << std::endl
      << " - middle box: "
      <<   (box_length - floor((mbox_max[0]-mbox_min[0])/dx_m)*dx_m) 
      <<   ", dx = " << dx_m << ", mismatch/dx = " 
      <<   (box_length/dx_m - floor((mbox_max[0]-mbox_min[0])/dx_m)) 
      << std::endl
      << " - bottom box: "
      <<   (box_length - floor((bbox_max[0]-bbox_min[0])/dx_t)*dx_t) 
      <<   ", dx = " << dx_t << ", mismatch/dx = " 
      <<   (box_length/dx_t - floor((bbox_max[0]-bbox_min[0])/dx_t)) 
      << std::endl;

  if constexpr (gdimension >= 3) 
    clog_one(warn) << "Lattice mismatch, Z-direction:" << std::endl
      << " -    top box: "
      <<   (box_height - floor((tbox_max[2]-tbox_min[2])/dz_t)*dz_t) 
      <<   ", dz = " << dz_t << ", mismatch/dz = " 
      <<   (box_height/dz_t - floor((bbox_max[2]-bbox_min[2])/dz_t)) 
      << std::endl
      << " - middle box: "
      <<   (box_height - floor((mbox_max[2]-mbox_min[2])/dz_m)*dz_m) 
      <<   ", dz = " << dz_m << ", mismatch/dz = " 
      <<   (box_height/dz_m - floor((mbox_max[2]-mbox_min[2])/dz_m)) 
      << std::endl
      << " - bottom box: "
      <<   (box_height - floor((bbox_max[2]-bbox_min[2])/dz_t)*dz_t) 
      <<   ", dz = " << dz_t << ", mismatch/dz = " 
      <<   (box_height/dz_t - floor((bbox_max[2]-bbox_min[2])/dz_t)) 
      << std::endl;

  // count the number of particles
  np_middle = particle_lattice::count(lattice_type,gdimension,mbox_min,mbox_max,
                                      sph_separation,0);
  np_top    = particle_lattice::count(lattice_type,gdimension,tbox_min,tbox_max,
                                      sph_sep_t, np_middle);
  np_bottom = particle_lattice::count(lattice_type,gdimension,bbox_min,bbox_max,
                                      sph_sep_t, np_middle + np_top);

  SET_PARAM(nparticles, np_middle + np_bottom + np_top);

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace param;

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

  // anything other than 2D is not implemented yet
  assert (gdimension == 2 || gdimension == 3);
  assert (domain_type == 0);

  // screen output
  clog_one(info)
    << "Kelvin-Helmholtz instability initial data " 
    << "in " << gdimension << "D" << std::endl;

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // screen output
  clog_one(info)
    << "Number of particles: "
    << nparticles << std::endl;
  clog_one(info)
    << "Initial data file: " 
    << initial_data_file << std::endl;

  // allocate arrays
  // Position
  double* x = new double[nparticles]();
  double* y = new double[nparticles]();
  double* z = new double[nparticles]();
  // Velocity
  double* vx = new double[nparticles]();
  double* vy = new double[nparticles]();
  double* vz = new double[nparticles]();
  // Acceleration
  double* ax = new double[nparticles]();
  double* ay = new double[nparticles]();
  double* az = new double[nparticles]();
  // Smoothing length
  double* h = new double[nparticles]();
  // Density
  double* rho = new double[nparticles]();
  // Internal Energy
  double* u = new double[nparticles]();
  // Pressure
  double* P = new double[nparticles]();
  // Mass
  double* m = new double[nparticles]();
  // Id
  int64_t* id = new int64_t[nparticles]();
  // Timestep
  double* dt = new double[nparticles]();

  // generate the lattice
  assert (np_middle == particle_lattice::generate( lattice_type,gdimension,
          mbox_min,mbox_max,sph_separation,0,x,y,z));
  assert (np_top    == particle_lattice::generate( lattice_type,gdimension,
          tbox_min,tbox_max,sph_sep_t,np_middle,x,y,z));
  assert (np_bottom == particle_lattice::generate( lattice_type,gdimension,
          bbox_min,bbox_max,sph_sep_t,nparticles-np_bottom,x,y,z));

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*std::max(pressure_m/rho_m,pressure_t/rho_t));

  // suggested timestep
  double timestep = timestep_cfl_factor
                  * sph_separation/std::max(cs, flow_velocity);

  // particle id number
  int64_t posid = 0;
  for(int64_t part=0; part<nparticles; ++part){
    id[part] = posid++;
    if (particle_lattice::in_domain_1d(y[part],
        mbox_min[1], mbox_max[1], domain_type)) {
      P[part] = pressure_m;
      rho[part] = rho_m;
      vx[part] = vx_m;
      m[part] = pmass;
    }
    else {
      P[part] = pressure_t;
      rho[part] = rho_t;
      vx[part] = vx_t;
      m[part] = pmass_t;
    }

    vy[part] = 0.;

    // Add velocity perturbation a-la Price (2008)
    if(y[part] < 0.25 and y[part] > 0.25 - 0.025)
      vy[part] = KH_A*sin(-2*M_PI*(x[part]+.5)/KH_lambda);
    if(y[part] >-0.25 and y[part] <-0.25 + 0.025)
      vy[part] = KH_A*sin(2*M_PI*(x[part]+.5)/KH_lambda);

    // compute internal energy using gamma-law eos
    u[part] = P[part]/(poly_gamma-1.)/rho[part];

    // particle masses and smoothing length
    m[part] = pmass;
    h[part] = sph_eta * kernels::kernel_width
                      * pow(m[part]/rho[part],1./gdimension);

  } // for part=0..nparticles

  // delete the output file if exists
  remove(initial_data_file.c_str());

  hid_t dataFile = H5P_openFile(initial_data_file.c_str(),H5F_ACC_RDWR);

  int use_fixed_timestep = 1;
  // add the global attributes
  H5P_writeAttribute(dataFile,"nparticles",&nparticles);
  H5P_writeAttribute(dataFile,"timestep",&timestep);
  int dim = gdimension;
  H5P_writeAttribute(dataFile,"dimension",&dim);
  H5P_writeAttribute(dataFile,"use_fixed_timestep",&use_fixed_timestep);

  H5P_setNumParticles(nparticles);
  H5P_setStep(dataFile,0);

  //H5PartSetNumParticles(dataFile,nparticles);
  H5P_writeDataset(dataFile,"x",x,nparticles);
  H5P_writeDataset(dataFile,"y",y,nparticles);
  H5P_writeDataset(dataFile,"z",z,nparticles);
  H5P_writeDataset(dataFile,"vx",vx,nparticles);
  H5P_writeDataset(dataFile,"vy",vy,nparticles);
  H5P_writeDataset(dataFile,"h",h,nparticles);
  H5P_writeDataset(dataFile,"rho",rho,nparticles);
  H5P_writeDataset(dataFile,"u",u,nparticles);
  H5P_writeDataset(dataFile,"P",P,nparticles);
  H5P_writeDataset(dataFile,"m",m,nparticles);
  H5P_writeDataset(dataFile,"id",id,nparticles);

  H5P_closeFile(dataFile);

  delete[] x, y, z, vx, vy, vz, ax, ay, az, h, rho, u, P, m, id, dt;

  MPI_Finalize();
  return 0;
}
