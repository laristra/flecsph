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
      << "Initial data generator for Rayleigh-Taylor (RT) instability test "
      << "in " << gdimension << "D" << std::endl
      << "Usage: ./RT_Xd_generator <parameter-file.par>" << std::endl;
}

//
// derived parameters
//
static double pressure_0;             // Initial pressure
static int64_t nparticlesproc;        // number of particles per proc
static double rho_1, rho_2;           // densities
static double vx_1, vx_2;             // velocities
static double pressure_1, pressure_2; // pressures
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

// geometric extents of the two regions: top and bottom
static point_t tbox_min, tbox_max;
static point_t bbox_min, bbox_max;

static int64_t np_top = 0;    // number of particles in the top block
static int64_t np_bottom = 0; // number of particles in the bottom block
static double sph_sep_t = 0; // particle separation in top or bottom blocks
static double pmass = 0;      // particle mass in the middle block

double u_from_eos(const double& rho, const double& p) {
  return p / ((param::poly_gamma-1.0)*rho);
}

void set_derived_params() {
  using namespace std;
  using namespace param;
  const double b_tol = particle_lattice::b_tol;

  bbox_max[0] = tbox_max[0] = box_length/2.;
  bbox_min[0] = tbox_min[0] =-box_length/2.;

  bbox_min[1] =  -box_width/2.;
  bbox_max[1] =  0.;
  tbox_min[1] =  0.;
  tbox_max[1] =  box_width/2.;

  if(gdimension == 3){
    bbox_max[2] = tbox_max[2] = box_height/2.;
    bbox_min[2] = tbox_min[2] =-box_height/2.;
  }

  pressure_0 = 2.5;

  // 1 = bottom 2 = top
  rho_1 = rho_initial;  // 2.0 by default
  rho_2 = rho_1 * density_ratio;

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

  // select particle lattice and kernel function
  particle_lattice::select();
  kernels::select();

  // particle mass and spacing
  SET_PARAM(sph_separation, (box_length/(double)(lattice_nx-1)));
  if(gdimension == 3){
    pmass = rho_1*sph_separation*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= 1./sqrt(2.);
    sph_sep_t = sph_separation*cbrt(rho_1/rho_2);
  }
  if(gdimension == 2){
    pmass = rho_1*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= sqrt(3)/2;
    sph_sep_t = sph_separation*sqrt(rho_1/rho_2);
  }

  // lattice spacing
  double dx, dy, dz, dx_t, dy_t, dz_t;
  dx = dy = dz = sph_separation;
  dx_t = dy_t = dz_t = sph_sep_t;
  if (lattice_type == 0) {
    clog_one(info)
      << "Lattice: rectangular, resolution: " << std::endl
      << " - top box:    dx = " << dx_t << std::endl
      << " - bottom box: dx = " << dx   << std::endl;
  }
  if (lattice_type == 1) { // HCP lattice
    dy   *= sqrt(3.);
    dy_t *= sqrt(3.);
    dz   *= 2.*sqrt(2./3.);
    dz_t *= 2.*sqrt(2./3.);
    clog_one(info)
      << "Lattice: HCP, resolution: " << std::endl
      << " - top box:     dx = " << dx_t << std::endl
      << "              2*dy = " << dy_t << std::endl
      << "              2*dz = " << dz_t << std::endl
      << " - bottom box:  dx = " << dx << std::endl
      << "              2*dy = " << dy << std::endl
      << "              2*dz = " << dz << std::endl;
  }
  if (lattice_type == 2) { // FCC lattice
    dy   *= sqrt(3.);
    dy_t *= sqrt(3.);
    dz   *= 3.*sqrt(2./3.);
    dz_t *= 3.*sqrt(2./3.);
    clog_one(info)
      << "Lattice: FCC, resolution: " << std::endl
      << " - top box:     dx = " << dx_t << std::endl
      << "              2*dy = " << dy_t << std::endl
      << "              3*dz = " << dz_t << std::endl
      << " - bottom box:  dx = " << dx << std::endl
      << "              2*dy = " << dy << std::endl
      << "              3*dz = " << dz << std::endl;
  }

  // adjust bottom lattice block in vertical direction
  bbox_min[1] += (box_width/2.) 
         - ((int)(box_width/2./dy))*dy;

  // for periodic boundaries, lattice has to match up:
  // adjust domain length
  if (periodic_boundary_x) {
    int Nx = (int)(box_length/dx) - 1;
    for(int i = Nx; i < Nx*100; ++i) {
      double w2 = floor(i*dx/dx_t + b_tol)*dx_t;
      if (fabs(w2 - i*dx) < lattice_mismatch_tolerance*dx_t) {
        SET_PARAM(box_length, std::min(w2,i*dx));
        bbox_min[0] = -box_length/2.;
        bbox_max[0] =  box_length/2.;
        tbox_min[0] = -box_length/2.;
        tbox_max[0] =  box_length/2.;
        break;
      }
    }
  }
  
  // adjust domain height
  if (gdimension >= 3 and periodic_boundary_z) {
    int Nz = (int)(box_length/dz) - 1;
    for(int k = Nz; k < Nz*100; ++k) {
      double w2 = floor(k*dz/dz_t + b_tol)*dz_t;
      if (fabs(w2 - k*dz) < lattice_mismatch_tolerance*dz_t) {
        SET_PARAM(box_height, std::min(w2,k*dz));
        bbox_min[2] = -box_height/2.;
        bbox_max[2] =  box_height/2.;
        tbox_min[2] = -box_height/2.;
        tbox_max[2] =  box_height/2.;
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

    clog_one(warn) << "Lattice mismatch, X-direction:" << std::endl
        << " -    top box: "
        <<   (box_length - floor((tbox_max[0]-tbox_min[0])/dx_t)*dx_t) 
        <<   ", dx = " << dx_t << ", mismatch/dx = " 
        <<   (box_length/dx_t - floor((tbox_max[0]-tbox_min[0])/dx_t)) 
        << std::endl
        << " - bottom box: "
        <<   (box_length - floor((bbox_max[0]-bbox_min[0])/dx)*dx) 
        <<   ", dx = " << dx << ", mismatch/dx = " 
        <<   (box_length/dx - floor((bbox_max[0]-bbox_min[0])/dx)) 
        << std::endl;

    if constexpr (gdimension >= 3) 
      clog_one(warn) << "Lattice mismatch, Z-direction:" << std::endl
        << " -    top box: "
        <<   (box_height - floor((tbox_max[0]-tbox_min[0])/dz_t)*dz_t) 
        <<   ", dz = " << dz_t << ", mismatch/dz = " 
        <<   (box_height/dz_t - floor((tbox_max[0]-tbox_min[0])/dz_t)) 
        << std::endl
        << " - bottom box: "
        <<   (box_height - floor((bbox_max[0]-bbox_min[0])/dz)*dz) 
        <<   ", dz = " << dz << ", mismatch/dz = " 
        <<   (box_height/dz - floor((bbox_max[0]-bbox_min[0])/dz)) 
        << std::endl;
  }

  // count the number of particles
  np_bottom = particle_lattice::count(lattice_type,gdimension,bbox_min,bbox_max,
                                      sph_separation, 0);
  np_top    = particle_lattice::count(lattice_type,gdimension,tbox_min,tbox_max,
                                      sph_sep_t, np_bottom);


  SET_PARAM(nparticles, np_bottom + np_top);

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
  //assert (gdimension == 2);
  assert (domain_type == 0);

  // screen output
  clog_one(info)
    << "Rayleigh-Taylor instability initial data " 
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
  assert (np_bottom == particle_lattice::generate( lattice_type,gdimension,
          bbox_min,bbox_max,sph_separation,0,x,y,z));
  assert (np_top    == particle_lattice::generate( lattice_type,gdimension,
          tbox_min,tbox_max,sph_sep_t,np_bottom,x,y,z));

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*std::max(pressure_1/rho_1,pressure_2/rho_2));

  // suggested timestep
  double timestep = timestep_cfl_factor
                  * sph_separation/std::max(cs, flow_velocity);

  // particle id number
  int64_t posid = 0;
  for(int64_t part=0; part<nparticles; ++part){
    id[part] = posid++;
    if (particle_lattice::in_domain_1d(y[part],
        bbox_min[1], bbox_max[1], domain_type)) {
      rho[part] = rho_1;
      m[part] = pmass;
      P[part] = pressure_0 
              + gravity_acceleration_constant
              *(rho_2*(tbox_max[1] - tbox_min[1]) - rho_1*y[part]);
    }
    else {
      rho[part] = rho_2;
      m[part] = pmass;
      P[part] = pressure_0 
              + gravity_acceleration_constant
              * rho_2*(tbox_max[1] - y[part]);
    }
    u[part] = u_from_eos(rho[part],P[part]);

    vx[part] = 0.;
    vy[part] = 0.;

    // Add velocity perturbation a-la Price (2008)
    if(fabs(y[part]) < .5*rt_perturbation_stripe_width) {
      if constexpr (gdimension == 2)
        vy[part] =-rt_perturbation_amplitude
                  *(1 + cos(2*M_PI*x[part]/box_length*rt_perturbation_mode))
                  *cos(M_PI*y[part]/rt_perturbation_stripe_width);

      if constexpr (gdimension == 3)
        vy[part] =-rt_perturbation_amplitude
                  *(1 + cos(2*M_PI*x[part]/box_length*rt_perturbation_mode))
                  *(1 + cos(2*M_PI*z[part]/box_length*rt_perturbation_mode))
                  *cos(M_PI*y[part]/rt_perturbation_stripe_width);
    }

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
