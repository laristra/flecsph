/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>
//#include <H5hut.h>

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
      << "Usage: ./RT_XD_generator <parameter-file.par>" << std::endl;
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

double pressure_gravity(const double& y, const double& rho) {
  using namespace param;
  return pressure_0 - rho * gravity_acceleration_constant * y;
}

double u_from_eos(const double& rho, const double& p) {
  return p / ((param::poly_gamma-1.0)*rho);
}

void set_derived_params() {
  using namespace std;
  using namespace param;

  // boundary tolerance factor
  const double b_tol = 1e-8;

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
  rho_2 = rho_1 * KH_density_ratio;

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

  std::cout<<"Boxes: " << std::endl << "up="
    <<tbox_min<<"-"<<tbox_max<<std::endl
    <<" bottom="
    <<bbox_min<<"-"<<bbox_max<< std::endl;

  // select particle lattice and kernel function
  particle_lattice::select();
  kernels::select();

  // particle mass and spacing
  SET_PARAM(sph_separation, (box_length*(1.0-b_tol)/(double)(lattice_nx-1)));
  if(gdimension == 3){
    pmass = rho_1*sph_separation*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= 1./sqrt(2.);
  }
  if(gdimension == 2){
    pmass = rho_1*sph_separation*sph_separation;
    if (lattice_type == 1 or lattice_type == 2)
      pmass *= sqrt(3)/2;
  }

  // adjust width of the middle block for symmetry
  double dy = sph_separation;
  if (lattice_type == 1 or lattice_type == 2)
    dy *= sqrt(3)/2;
  double dy_tb = dy; // lattice step in y-direction for top and bottom blocks
  //double bbox_width = bbox_max[1] - bbox_min[1];
  //bbox_width = (int)(bbox_width/(2*dy))*2*dy;
  //bbox_min[1] = -bbox_width/2.;
  //bbox_max[1] =  bbox_width/2.;

  sph_sep_t = sph_separation * sqrt(rho_1/rho_2);

  //dy_tb = dy * sph_sep_t/sph_separation;

  // adjust top blocks
  tbox_min[1] = bbox_max[1] - dy + 0.5*(dy_tb + dy);

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

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // screen output
  std::cout << "Kelvin-Helmholtz instability setup in " << gdimension
       << "D:" << std::endl << " - number of particles: " << nparticles << std::endl
       << " - generated initial data file: " << initial_data_file << std::endl;

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
    }
    else {
      rho[part] = rho_2;
      m[part] = pmass;
    }

    P[part] = pressure_gravity(y[part],rho[part]);
    u[part] = u_from_eos(rho[part],P[part]);

    vy[part] = 0.;

    // Add velocity perturbation a-la Price (2008)
    //vy[part] = 0.01*(1 + cos(4*M_PI*x[part]))*(1 + cos(3*M_PI*y[part]))/4.;
    if(y[part] < 0.025 and y[part] > -0.025)
      vy[part] = 0.01*cos(M_PI*(x[part]/box_length));

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
