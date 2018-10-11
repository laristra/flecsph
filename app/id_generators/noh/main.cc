/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


/*
 * Noh Collapse test
 * -----------------
 * The test is initialized as a disk / sphere of particles with homogeneous 
 * density, vanishingly small pressure, and inward velocity v_r. As particles 
 * move inwards, they pile up at the center at forming a region with constant 
 * density that is dependent on gamma and the dimensionality of the problem.
 * A standing shock front forms that moves outwards as more particles are piling
 * up. Its radial distance grows as 
 *
 * r_shock (t) = 0.5 (gamma - 1) v_r t
 *
 * The density of infalling matter evolves as:
 *
 * n(r>=r_shock) = n_0 (1 + (v_r/r)t)^(d-1)
 *
 * The density of matter enclose by the shock is given by:
 *
 * n(r<r_rhock) = n_0 * ((gamma + 1)/(gamma - 1))^d
 *
 * where d gives the geometry of the system (1=planar, 2=cyl., 3=spherical)
 *
 * For more information, see:
 * Liska & Wendroff, 2003, SIAM J. Sci. Comput. 25, 995, Section 4.5
*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>
#include <H5hut.h>

#include "user.h"
#include "noh.h"
#include "params.h"

//
// help message
//
void print_usage() {
  clog(warn)
      << "Initial data generator for the Noh collapse test" << std::endl
      << "Usage: ./noh_generator <parameter-file.par>"      << std::endl;
}


//
// derived parameters
//
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  // total number of particles
  SET_PARAM(nparticles, lattice_nx*lattice_nx);
  SET_PARAM(sph_smoothing_length, (2.0*sph_separation));

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}


// Is particle with x,y coordinates within a disk (center located at x0,y0) with radius r?
bool in_radius(
    double x,
    double y,
    double x0,
    double y0,
    double r) {
    return (x-x0)*(x-x0)+(y-y0)*(y-y0) < r*r;
}


int main(int argc, char * argv[]){
  using namespace param;

  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // check options list: exactly one command-line argument is accepted
  if (argc != 2) {
    clog(error) << "ERROR: parameter file not specified!" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // Radius of the disk given the number of particles
  double radius = (sph_separation*(lattice_nx - 1))/2.0;

  // Center of disk
  double x_c = (lattice_nx - 1)*sph_separation/2.0;
  double y_c = x_c;


  clog(info) << "Sphere: r=" << radius << std::endl
                 << "origin: pos=["<<x_c<<";"<<y_c<<"]" << std::endl
                 << "Generating "  << nparticles
                 << " particles (="<< lattice_nx<<"^2) " << std::endl;


  // Start on  0 0
  double x_topproc = x_c - radius;
  double y_topproc = y_c - radius;

  double maxxposition = x_c + radius;
  double maxyposition = y_c + radius;

  // Particle positions
  double* x = new double[nparticles]();
  double* y = new double[nparticles]();
  double* z = new double[nparticles]();

  // Particle velocities
  double* vx = new double[nparticles]();
  double* vy = new double[nparticles]();
  double* vz = new double[nparticles]();

  // Particle accelerations
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

  // ID
  int64_t* id = new int64_t[nparticles]();

  // Timestep
  double* dt = new double[nparticles]();

  // ID of first particle
  int64_t posid = 0;


  // Header data
  // the number of particles = nparticles
  // The value for constant timestep
  double timestep = initial_dt;  // TODO: replace with parameter initial_dt
  int dimension = 2;             // TODO: use global dimension parameter

  clog(info) << "top_X=" << x_topproc << " top_Y="  << y_topproc    << std::endl
                 << "maxX=" << maxxposition << " maxY=" << maxyposition << std::endl;

  double xposition = x_topproc;
  double yposition = y_topproc;
  int64_t tparticles = 0;
  double max_radius = 0;

  // Loop over all particles and assign position, velocity etc.
  for (int64_t part = 0; part < nparticles; ++part) {

    // Checks if particle is in the disk and creates position coordinates
    while (!in_radius(xposition, yposition, x_c, y_c, radius)) {
      xposition += sph_separation;
      if (xposition > maxxposition) {
        if (yposition > maxyposition) {break;}
        xposition  = x_topproc;
        yposition += sph_separation;
      }
    }
    if (xposition > maxxposition) {
      if (yposition > maxyposition) {break;}
    }

    tparticles++;

    // Assign particle position
    x[part] = xposition;
    y[part] = yposition;

    // Determine the maximum radial distance for a particle for 
    // later mass/density determination 
    if (sqrt((x[part]-x_c)*(x[part]-x_c) + (y[part]-y_c)*(y[part]-y_c)) 
    > max_radius) {
        max_radius = sqrt((x[part]-x_c)*(x[part]-x_c) 
                        + (y[part]-y_c)*(y[part]-y_c));
    }

    xposition += sph_separation;
    if (xposition > maxxposition) {
      if (yposition > maxyposition) {break;}
      xposition  = x_topproc;
      yposition += sph_separation;
    }

    // Assign particle pressure
    P[part] = pressure_initial;

    // Assign particle density
    rho[part] = rho_initial;

    // Assign particle accelerations
    ax[part] = 0.0;
    ay[part] = 0.0;
    az[part] = 0.0;

    // Assign particle mass
    //m[part] = m_in;

    // Assign particle internal energy
    u[part] = P[part]/(poly_gamma-1.)/rho[part];

    // Assign particle smoothing length
    h[part] = sph_smoothing_length;

    // Assign particle ID
    id[part] = posid++;

    // Assign particle inward pointing velocity with absolute value 0.1
    double A = sqrt((x[part]-x_c)*(x[part]-x_c) + (y[part]-y_c)*(y[part]-y_c));
    if (A <= 0.0) {
      vx[part] = 0.0;
      vy[part] = 0.0;
    }
    else {
      vx[part] = -(x[part]-x_c) * 0.1 / A;
      vy[part] = -(y[part]-y_c) * 0.1 / A;
    }
  }

  clog(info) << "Actual number of particles inside the sphere: "
                 << tparticles << std::endl;

  // Particle mass given the density, number of particles, and disk radius
  double m_in = rho_initial * max_radius * max_radius * M_PI / tparticles;

  for (int64_t part = 0; part < nparticles; ++part) m[part] = m_in;


  // remove the previous file
  remove(initial_data_file.c_str());

  h5_file_t * dataFile = H5OpenFile(initial_data_file.c_str()
      ,H5_O_WRONLY, MPI_COMM_WORLD);
    
  int use_fixed_timestep = 1; 
  // add the global attributes
  H5WriteFileAttribInt64(dataFile,"nparticles",&tparticles,1);
  H5WriteFileAttribFloat64(dataFile,"timestep",&timestep,1);
  int dim = gdimension;
  H5WriteFileAttribInt32(dataFile,"dimension",&dim,1);
  H5WriteFileAttribInt32(dataFile,"use_fixed_timestep",&use_fixed_timestep,1);

  H5SetStep(dataFile,0);
  H5PartSetNumParticles(dataFile,tparticles);
  H5PartWriteDataFloat64(dataFile,"x",x);
  H5PartWriteDataFloat64(dataFile,"y",y);
  H5PartWriteDataFloat64(dataFile,"z",z);
  H5PartWriteDataFloat64(dataFile,"vx",vx);
  H5PartWriteDataFloat64(dataFile,"vy",vy);
  H5PartWriteDataFloat64(dataFile,"vz",vz);
  H5PartWriteDataFloat64(dataFile,"ax",ax);
  H5PartWriteDataFloat64(dataFile,"ay",ay);
  H5PartWriteDataFloat64(dataFile,"az",az);
  H5PartWriteDataFloat64(dataFile,"h",h);
  H5PartWriteDataFloat64(dataFile,"rho",rho);
  H5PartWriteDataFloat64(dataFile,"u",u);
  H5PartWriteDataFloat64(dataFile,"P",P);
  H5PartWriteDataFloat64(dataFile,"m",m);
  H5PartWriteDataInt64(dataFile,"id",id);
 
  H5CloseFile(dataFile);

  delete[] x, y, z, vx, vy, vz, ax, ay, az, h, rho, u, P, m, id, dt; 

  MPI_Finalize();
  return 0;
}
