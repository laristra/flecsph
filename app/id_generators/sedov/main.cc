/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>
#include <random>
#include <math.h>

#include "user.h"
#include "sedov.h"
#include "params.h"
#include "lattice.h"
#include "kernels.h"
#include "density_profiles.h"
#include "io.h"
using namespace io;
#include "bodies_system.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

/*
The Sedov test is set up with uniform density and vanishingly small pressure.
An explosion is initialized via a point-like deposition of energy E_blast in
the center of the simulation space. The resulting spherically symmetric shock
waves moves outwards with a radial distance given by

r_shock = (E_blast * t^2 / (alpha * n_0))^(1/(2*d))

where alpha is a constant of the order one with its exact value given by the
adiabatic index gamma. The peak density of the shock is given by:

n_shock = n_0 ((gamma + 1)/(gamma -1))

The velocity of shocked matter has a radial dependence approximately ~ r/t.
As the shock wave moves away from the center it leaves behind matter at
vanishingly low density. With the pressure staying finite for r=0, the
temperature grows and becomes infinitely large at the origin of the blast wave .

Reference:
G.  Taylor, “The Formation of a Blast Wave by a Very Intense Explosion.
I.  Theoretical  Discussion,” Royal Society of London Proceedings Series A
, vol. 201, pp. 159–174, Mar. 1950.
*/


//
// help message
//
void print_usage() {
  std::cout
      << "Initial data generator for the " << gdimension << "D Sedov blast wave"
      << std::endl << "Usage: ./sedov_generator <parameter-file.par>"
      << std::endl;
}

//
// derived parameters
//
static double timestep = 1.0;         // Recommended timestep
static double total_mass = 1.;        // total mass of the fluid
static double mass_particle = 1.;     // mass of an individual particle
static point_t bbox_max, bbox_min;    // bounding box of the domain
static char initial_data_file[256];   // = initial_data_prefix[_XXXXX].h5part"

void set_derived_params() {
  using namespace param;

  density_profiles::select();
  particle_lattice::select();

  // The value for constant timestep
  timestep = initial_dt;

  // Bounding box of the domain
  if (domain_type == 0) { // box
    bbox_min[0] = -box_length/2.;
    bbox_max[0] =  box_length/2.;
    if constexpr (gdimension > 1) {
      bbox_min[1] = -box_width/2.;
      bbox_max[1] =  box_width/2.;
    }
    if constexpr (gdimension > 2) {
      bbox_min[2] = -box_height/2.;
      bbox_max[2] =  box_height/2.;
    }
  }
  else if (domain_type == 1) { // sphere or circle
    bbox_min = -sphere_radius;
    bbox_max =  sphere_radius;
  }

  // particle separation
  if (domain_type == 0) {
    SET_PARAM(sph_separation, (box_length/(lattice_nx-1)));
  }
  else if (domain_type == 1) {
    SET_PARAM(sph_separation, (2.*sphere_radius/(lattice_nx-1)));
  }

  // Count number of particles
  int64_t tparticles =
      particle_lattice::count(lattice_type,domain_type,
      bbox_min,bbox_max,sph_separation,0);
  SET_PARAM(nparticles, tparticles);

  // total mass
  if constexpr (gdimension == 1) {
    total_mass = rho_initial * box_length;
  }
  if constexpr (gdimension == 2) {
    if (domain_type == 0) // a box
      total_mass = rho_initial * box_length*box_width;
    else if (domain_type == 1) // a circle
      total_mass = rho_initial * M_PI*SQ(sphere_radius);
  }
  if constexpr (gdimension == 3) {
    if (domain_type == 0) { // a box
      assert (boost::iequals(density_profile,"constant"));
      total_mass = rho_initial * box_length*box_width*box_height;
    }
    else if (domain_type == 1) { // a sphere
      // normalize mass such that central density is rho_initial
      total_mass = rho_initial * CU(sphere_radius)
                 / density_profiles::spherical_density_profile(0.0);
    }
  }

  // single particle mass
  assert (equal_mass);
  mass_particle = total_mass / nparticles;

  // set kernel
  kernels::select();

  // smoothing length
  const double sph_h = sph_eta * kernels::kernel_width
                               * pow(mass_particle/rho_initial,1./gdimension);
  SET_PARAM(sph_smoothing_length, sph_h);

  // intial internal energy
  SET_PARAM(uint_initial, (pressure_initial/(rho_initial*(poly_gamma-1.0))));

  // Filename to be generated
  bool input_single_file = H5P_fileExists(initial_data_prefix);  
  if (input_single_file or initial_iteration == 0) 
    sprintf(initial_data_file,"%s.h5part",initial_data_prefix);
  else {
    
    // find the file with initial_iteration
    int step = H5P_findIterationSnapshot(initial_data_prefix, 
                                  param::initial_iteration);
    // file doesn't exist: complain and exit
    if (step < 0) {
      clog(error) << "Cannot find iteration " << param::initial_iteration 
                  <<" in prefix " << initial_data_prefix << std::endl;
      exit(MPI_Barrier(MPI_COMM_WORLD) && MPI_Finalize());
    }
    sprintf(initial_data_file,"%s_%05d.h5part",initial_data_prefix,step);
  }

}

int main(int argc, char * argv[]){
  using namespace param;

  // check options list: exactly one option is allowed
  if (argc != 2) {
    std::cerr << "ERROR: parameter file not specified!" << std::endl;
    print_usage();
    exit(0);
  }

  // launch MPI
  int rank, size;
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  assert (size == 1); // parallel ID generator not implemented yet
  clog_set_output_rank(0);

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();
  body_system<double,gdimension> bs;
  if (modify_initial_data) { 
    bs.read_bodies(initial_data_prefix,"",initial_iteration);
    SET_PARAM(nparticles, bs.getNBodies());
  }
  else {
    bs.getLocalbodies().clear();
    bs.getLocalbodies().resize(nparticles);
  }
  auto & bodies = bs.getLocalbodies();

  // Declare coordinate arrays
  double* x = new double[nparticles]();
  double* y = new double[nparticles]();
  double* z = new double[nparticles]();

  if (not modify_initial_data) { 
    // Generate the lattice
    assert(nparticles ==
        particle_lattice::generate(lattice_type,domain_type,
        bbox_min,bbox_max,sph_separation,0, x, y, z));

    for (int64_t a=0L; a<nparticles; ++a) {
      body& particle = bodies[a];
      if constexpr (gdimension == 1) {
        point_t pos = {x[a]};
        particle.set_coordinates(pos);
      }

      if constexpr (gdimension == 2) {
        point_t pos = {x[a],y[a]};
        particle.set_coordinates(pos);
      }

      if constexpr (gdimension == 3) {
        point_t pos = {x[a],y[a],z[a]};
        particle.set_coordinates(pos);
      }
    }
  }

  // Number of particles in the blast zone
  int64_t particles_blast = 0;

  // Total mass of particles in the blast zone
  double mass_blast = 0;

  // Count the number of particles and mass in the blast zone
  // The blast is centered at the origin ({0,0} or {0,0,0})
  for(int64_t a=0L; a<nparticles; ++a) {
    body& particle = bodies[a];
    double r = norm2(particle.coordinates());
    if (r < sedov_blast_radius) {
       particles_blast++;
       mass_blast += mass_particle;
    }
  }

  // Assign density, pressure and specific internal energy to particles,
  // including the particles in the blast zone
  const double rho0 = density_profiles::spherical_density_profile(0);
  const double K0 = pressure_initial // polytropic constant
                  / pow(rho_initial, poly_gamma); 
  std::default_random_engine generator;
  for(int64_t a=0; a<nparticles; ++a){
    body& particle = bodies[a];

    // zero velocity for this test
    point_t zero = 0;
    particle.setVelocity(zero);
    particle.setAcceleration(zero);

    // radial distance from the origin
    point_t rp(particle.coordinates());
    double r = norm2(rp);

    // set density, particle mass, smoothing length and id
    double rho_a, m_a, h_a;
    int64_t id_a;
    if (modify_initial_data) {
      rho_a = particle.getDensity();
      m_a   = particle.mass();
      h_a   = particle.radius();
      id_a  = particle.id();
    }
    else {
      rho_a = rho_initial/rho0  // renormalize density profile
            * density_profiles::spherical_density_profile(r/sphere_radius);
      m_a = mass_particle;
      h_a = sph_eta * kernels::kernel_width
          * pow(mass_particle/rho_a,1./gdimension);
      id_a = a;
      particle.setDensity(rho_a);
      particle.set_mass(m_a);
      particle.set_radius(h_a);
      particle.set_id(a);
    }

    if (lattice_perturbation_amplitude > 0.0) {
    // add lattice perturbation
      std::normal_distribution<double> 
        distribution(0.,h_a*lattice_perturbation_amplitude);
      for (unsigned short k=0; k<gdimension; ++k) { 
        rp[k] += distribution(generator);
      }
      particle.set_coordinates(rp);
    }

    // Blast energy in input file is given as total energy.
    // FleCSPH uses specific internal energy. 
    // Convert blast energy in specific internal energy: 
    double u_blast = sedov_blast_energy/mass_blast;

    // set internal energy
    double u_a = K0*pow(rho_a,poly_gamma-1)/(poly_gamma-1);
    if (r < sedov_blast_radius) 
      u_a += u_blast;
      //u_a += sedov_blast_energy/particles_blast;
    particle.setInternalenergy(u_a);

    // set pressure (a function of density and internal energy)
    double P_a = rho_a*u_a*(poly_gamma - 1);
    particle.setPressure(P_a);

    // set timestep
    particle.setDt(initial_dt);
  }

  clog_one(info) << "Number of particles: " << nparticles << std::endl;
  clog_one(info) << "Total number of seeded blast particles: " << particles_blast << std::endl;
  clog_one(info) << "Mass of seeded blast particles: " << mass_blast << std::endl;
  clog_one(info) << "Total blast energy: " << sedov_blast_energy << std::endl;

  // remove the previous file
  remove(initial_data_file);
  delete[] x, y, z;

  // write the file; iteration for initial data MUST BE zero!!
  bs.write_bodies(initial_data_prefix, 0, 0.0);
  MPI_Finalize();
  return 0;
}
