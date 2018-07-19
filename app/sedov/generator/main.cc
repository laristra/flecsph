/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>

#include "params.h"
#include "hdf5ParticleIO.h"
#include "kernels.h"
#include "user.h"

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
  clog_one(warn)
      << "Initial data generator for the 2D Sedov blast wave" << std::endl
      << "Usage: ./sedov_generator <parameter-file.par>"      << std::endl;
}

//
// derived parameters
//
static double r_blast;      // Radius of injection region
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace param;

  // total number of particles
  SET_PARAM(nparticles, sqrt_nparticles*sqrt_nparticles);
  SET_PARAM(uint_initial, (pressure_initial/(rho_initial*(poly_gamma-1.0))));
  SET_PARAM(sph_smoothing_length, (5.*sph_separation));
  r_blast = sedov_blast_radius * sph_separation;   // Radius of injection region

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}


bool
in_radius(
    double x,
    double y,
    double z,
    double x0,
    double y0,
    double z0,
    double r)
{
  return (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0)<r*r;
}

int main(int argc, char * argv[]){
  using namespace param;

  // launch MPI
  int rank, size;
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // check options list: exactly one option is allowed
  if (argc != 2) {
    clog_one(error) << "ERROR: parameter file not specified!" << std::endl;
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // Header data
  // the number of particles = nparticles
  // The value for constant timestep
  double timestep = 0.001;  // TODO: replace with parameter initial_dt
  int64_t maximum_part = nparticles;
  if (gdimension==3){
    maximum_part = pow(sqrt_nparticles,3);
  }

  // Start on  0 0
  double radius = (sph_separation*(sqrt_nparticles-1.))/2.;

  // Central coordinates
  double x_c;
  double y_c;
  double z_c;
  if (gdimension == 2){
    x_c = (sqrt_nparticles-1)*sph_separation/2.;
    y_c = x_c;
    z_c = 0;
  } else if(gdimension == 3){
    x_c = (sqrt_nparticles-1)*sph_separation/2.;
    y_c = x_c;
    z_c = x_c;
  }

  // Maximum coordinates of the initial configuration
  double maxxposition;
  double maxyposition;
  double maxzposition;
  if (gdimension == 2){
    maxxposition = x_c + radius;
    maxyposition = y_c + radius;
    maxzposition = 0;
  } else if(gdimension == 3){
    maxxposition = x_c + radius;
    maxyposition = y_c + radius;
    maxzposition = z_c + radius;
  }

  if (gdimension == 2){
    clog_one(info) <<"Sphere: r="<<radius<<" pos=["<<x_c<<";"<<y_c<<"]"<<std::endl
          << "Attempting to generate " << maximum_part << " particles" << std::endl;
  } else if(gdimension == 3){
    clog_one(info) <<"Sphere: r="<<radius<<" pos=["<<x_c<<";"<<y_c<<";"<<z_c<<"]"<<std::endl
          << "Attempting to generate " << maximum_part << " particles" << std::endl;
  }

  // Setting the coordinates to start building the lattice
  double x_topproc;
  double y_topproc;
  double z_topproc;
  if (gdimension == 2){
    x_topproc = x_c - radius;
    y_topproc = y_c - radius;
    z_topproc = 0;
  } else if(gdimension == 3){
    x_topproc = x_c - radius;
    y_topproc = y_c - radius;
    z_topproc = z_c - radius;
  }

  // Position
  double* x = new double[maximum_part]();
  double* y = new double[maximum_part]();
  double* z = new double[maximum_part]();
  // Velocity
  double* vx = new double[maximum_part]();
  double* vy = new double[maximum_part]();
  double* vz = new double[maximum_part]();
  // Acceleration
  double* ax = new double[maximum_part]();
  double* ay = new double[maximum_part]();
  double* az = new double[maximum_part]();
  // Smoothing length
  double* h = new double[maximum_part]();
  // Density
  double* rho = new double[maximum_part]();
  // Internal Energy
  double* u = new double[maximum_part]();
  // Pressure
  double* P = new double[maximum_part]();
  // Mass
  double* m = new double[maximum_part]();
  // Id
  int64_t* id = new int64_t[maximum_part]();
  // Timestep
  double* dt = new double[maximum_part]();

  // Id of my first particle
  int64_t posid = 0;

  // Number of particles in the blast zone
  int64_t particles_blast = 0;

  // Total mass of particles in the blast zone
  double mass_blast = 0;

  double xposition = x_topproc;//0;
  double yposition = y_topproc;//0;
  double zposition = z_topproc;//0;

  int64_t tparticles = 0;

  // For square lattice initial configuration
  if(lattice_type == 0) {
    for (int64_t part=0; part<maximum_part; ++part) {
      while(!in_radius(xposition,yposition,zposition,x_c,y_c,z_c,radius)){
        xposition+= sph_separation;
        if(xposition > maxxposition){
          if(yposition > maxyposition){
            if(gdimension==3){
              if(zposition > maxzposition){
                break;
              }
              zposition+=sph_separation;
              xposition=x_topproc;
              yposition=y_topproc;
            } else if(gdimension==2){
              break;
            }
          }
          xposition=x_topproc;
          yposition+=sph_separation;
        }
      }

      if(xposition > maxxposition){
        if(yposition > maxyposition){
          if(gdimension==3){
            if(zposition > maxzposition){
              break;
            }
          } else if(gdimension==2){
            break;
          }
        }
      }

      tparticles++;
      x[part] = xposition;
      y[part] = yposition;
      z[part] = zposition;

      xposition+=sph_separation;
      if(xposition > maxxposition){
        if(yposition > maxyposition){
          if(gdimension==3){
            if(zposition > maxzposition){
              break;
            }
            zposition+=sph_separation;
            xposition=x_topproc;
            yposition=y_topproc;
          } else if(gdimension==2){
            break;
          }
        }
        xposition=x_topproc;
        yposition+=sph_separation;
      }

    }
  } else if (lattice_type == 1){ // For triangular lattice initial configuration
    int64_t row = 0;
    int64_t zrow = 0;
    for (int64_t part=0; part<maximum_part; ++part) {
      while(!in_radius(xposition,yposition,zposition,x_c,y_c,z_c,radius)){
        xposition+= sph_separation;
        if(xposition > maxxposition){
          if(yposition > maxyposition){
            if(gdimension==3){
              if(zposition > maxzposition){
                break;
              }
              zposition+=sqrt(2./3.)*sph_separation;
              zrow+=1;
              if (zrow % 2 == 1){
                xposition=x_topproc-sph_separation/2.;
                yposition=y_topproc-sph_separation*sqrt(3.)/6.;
                row=0;
              } else if (zrow % 2 == 0){
                xposition=x_topproc;
                yposition=y_topproc;
                row=0;
              }

            } else if(gdimension==2){
              break;
            }
          }
          yposition+=sqrt(3.)/2.*sph_separation;
          row+=1;
          if (row % 2 == 1){
            xposition=x_topproc-sph_separation/2.;
          } else if (row % 2 == 0){
            xposition=x_topproc;
          }
        }
      }

      if(xposition > maxxposition){
        if(yposition > maxyposition){
          if(gdimension==3){
            if(zposition > maxzposition){
              break;
            }
          } else if(gdimension==2){
            break;
          }
        }
      }

      tparticles++;
      x[part] = xposition;
      y[part] = yposition;
      z[part] = zposition;

      xposition+=sph_separation;
      if(xposition > maxxposition){
        if(yposition > maxyposition){
          if(gdimension==3){
            if(zposition > maxzposition){
              break;
            }
            zposition+=sqrt(2./3.)*sph_separation;
            zrow+=1;
            if (zrow % 2 == 1){
              xposition=x_topproc-sph_separation/2.;
              yposition=y_topproc-sph_separation*sqrt(3.)/6.;
              row=0;
            } else if (zrow%2 == 0){
              xposition=x_topproc;
              yposition=y_topproc;
              row=0;
            }
          } else if(gdimension==2){
            break;
          }
        }
        yposition+=sqrt(3.)/2.*sph_separation;
        row += 1;
        if (row % 2 == 1){
          xposition=x_topproc-sph_separation/2.;
        } else if (row % 2 == 0){
          xposition=x_topproc;
        }
      }

    }
  }
  double mass = rho_initial * M_PI*pow(radius,2.)/tparticles;
  // Assign density, pressure and specific internal energy to particles,
  // including the particles in the blast zone
  for(int64_t part=0; part<tparticles; ++part){
    // Particle mass from number of particles and density
    m[part] = mass;

    // Count particles in the blast zone and sum their masses
    if(sqrt((x[part]-x_c)*(x[part]-x_c)+(y[part]-y_c)*(y[part]-y_c)+(z[part]-z_c)*(z[part]-z_c)) < r_blast){
       particles_blast++;
       mass_blast += m[part];
    }
    P[part] = pressure_initial;
    rho[part] = rho_initial;
    u[part] = uint_initial;
    h[part] = sph_smoothing_length;
    id[part] = posid++;

    if(sqrt((x[part]-x_c)*(x[part]-x_c)+(y[part]-y_c)*(y[part]-y_c)+(z[part]-z_c)*(z[part]-z_c)) < r_blast){
       u[part] = u[part]+sedov_blast_energy/particles_blast;
       P[part] = u[part]*rho[part]*(poly_gamma - 1.0);
    }
  }

  clog_one(info) << "Real number of particles: " << tparticles << std::endl;
  clog_one(info) << "Total number of seeded blast particles: " << particles_blast << std::endl;
  clog_one(info) << "Total blast energy (E_blast = u_blast * total mass): "
                 << sedov_blast_energy * mass_blast << std::endl;

  // remove the previous file
  remove(initial_data_file.c_str());

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(initial_data_file.c_str(),MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",tparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",gdimension);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",tparticles,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",tparticles,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",tparticles,z);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",tparticles,vx);
  _d2.createVariable("vy",Flecsi_Sim_IO::point,"double",tparticles,vy);
  _d3.createVariable("vz",Flecsi_Sim_IO::point,"double",tparticles,vz);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("ax",Flecsi_Sim_IO::point,"double",tparticles,ax);
  _d2.createVariable("ay",Flecsi_Sim_IO::point,"double",tparticles,ay);
  _d3.createVariable("az",Flecsi_Sim_IO::point,"double",tparticles,az);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();


  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",tparticles,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",tparticles,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",tparticles,u);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",tparticles,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",tparticles,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",tparticles,id);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  testDataSet.closeFile();

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] ax;
  delete[] ay;
  delete[] az;
  delete[] h;
  delete[] rho;
  delete[] u;
  delete[] P;
  delete[] m;
  delete[] id;
  delete[] dt;

  MPI_Finalize();
  return 0;
}
