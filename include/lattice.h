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
 * @file lattice_type.h
 * @author Alexander Kaltenborn
 * @date July 2018
 * @brief function to choose lattice type and populate initial conditions
 *
 * Module for generating lattices of different types, suitable for a wide variety
 * of problems, making it essential building block for initial data construction.
 *
 * Having dimension set in the user.h file, this function will currently populate
 * either a rectangular lattice, a triangular lattice, an FCC lattice, or a HCP
 * lattice.
 *
 * Inputs:
 *  lattice_type - int value, corresponding to 0:rect, 1:HCP, 2:FCC (1 & 2 are
 *                 triangular in 2D)
 *  domain_type  - int value, 0:cube, 1:sphere, 2:full domain
 *  bbox_min     - bottom corner location for smallest cube covering the domain
 *  bbox_max     - top corner location for smallest cube covering the domain
 *  sph_sep      - the separation for the given domain within bbox_min and bbox_max
 *  posid        - enter the particle ID number from which the function is called
 *                 it is important for when the function is assigning position
 *                 coordinates to the particles
 *  count_only   - boolean for determing particle number (returned for allocating
 *                 proper arrays) or writing positions to those arrays
 *  x,y,z[]      - the arrays to be filled for the positions of each particle
 *
 *
 * TODO: include variables for paralellization (mpi rank)
 */

#include <stdlib.h>
#include "user.h"
#include "tree.h"
#include <math.h>

namespace particle_lattice {

/**
 * @brief      in_domain checks to see if the entered particle position info
 *             is valid within the restrictive domain_type and total domain
 *             Returns boolean: true or false
 *
 * @param      x,y,z       - position of the particle
 *             x0,y0,z0    - position of the center of the domain (relevant for
 *                           spheres and cubes)
 *             bbox_min    - minimum position for the total domain
 *             bbox_max    - maximum position for the total domain
 *             r           - radius of the sphere or 1/2 length of cube edge
 *             domain_type - int value, 0:cube, 1:sphere, 2:full domain
 */
bool in_domain_1d(
    const double x,
    const double x0,
    const double xmin,
    const double xmax,
    const double r,
    const int domain_type)
{
  // within_domain checks to see if the position is within the total domain
  bool within_domain = (x>=xmin && x<=xmax);

  // Check if within domain_type
  if (domain_type == 0 || domain_type == 1)
    within_domain *= (std::abs(x-x0)<=r);
  return within_domain;
}

bool in_domain_2d(
    const double x,
    const double y,
    const double x0,
    const double y0,
    const point_t& bbox_min,
    const point_t& bbox_max,
    const double r,
    const int domain_type)
{
  // within_domain checks to see if the position is within the total domain
  double xmin = bbox_min[0], xmax = bbox_max[0];
  double ymin = bbox_min[1], ymax = bbox_max[1];
  bool within_domain = x>=xmin && x<=xmax && y>=ymin && y<=ymax;

  // Check if within domain_type
  if(domain_type==0)
    within_domain *= (std::abs(x-x0)<=r && std::abs(y-y0)<=r);
  else if(domain_type==1)
    within_domain *= ((x-x0)*(x-x0)+(y-y0)*(y-y0) <r*r);
  return within_domain;
}

bool in_domain_3d(
    const double x,
    const double y,
    const double z,
    const double x0,
    const double y0,
    const double z0,
    const point_t& bbox_min,
    const point_t& bbox_max,
    const double r,
    const int domain_type)
{
  // within_domain checks to see if the position is within the total domain
  double xmin = bbox_min[0], xmax = bbox_max[0];
  double ymin = bbox_min[1], ymax = bbox_max[1];
  double zmin = bbox_min[2], zmax = bbox_max[2];
  bool within_domain = 
       x>=xmin && x<=xmax && y>=ymin && y<=ymax && z>=zmin && z<=zmax;

  // Check if within domain_type
  if(domain_type==0)
    within_domain *= (std::abs(x-x0)<=r && std::abs(y-y0)<=r && std::abs(z-z0)<=r);
  else if(domain_type==1)
    within_domain *= ((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0) < r*r);
  return within_domain;
}

/**
 * @brief      Generate lattice will run through the supplied domain, and--
 *             depending on the count_only switch--count total number of particles
 *             or assign the positions to the position arrays
 *             Returns int64_t: total particle number
 *
 * @param      Refer to inputs section in introduction
 */
int64_t generator_lattice_1d(
   const int lattice_type,
   const int domain_type,
   const point_t& bbox_min,
   const point_t& bbox_max,
   const double sph_sep,
   int64_t posid,
   bool count_only = true,
   double x[] = NULL,
   double y[] = NULL,
   double z[] = NULL)
{
  // Determine radius of the domain as a radius of inscribed sphere
  double radius = 0.0;
  if(domain_type==1 || domain_type==0) {
     radius = 0.5*(bbox_max[0] - bbox_min[0]);
  }

  // Central coordinates: in most cases this should be centered at 0
  double x_c = (bbox_max[0] + bbox_min[0])/2.;

  // True number of particles to be determined and returned
  int64_t tparticles = 0;

  // regular lattice in 1D
  double xmin = bbox_min[0], xmax = bbox_max[0];
  for(double x_p = xmin; x_p < xmax; x_p += sph_sep){
    if(in_domain_1d(x_p,x_c,xmin,xmax,radius,domain_type)){
      tparticles++;
      if(!count_only){
        x[posid] = x_p;
        y[posid] = 0.0;
        z[posid] = 0.0;
      }
      posid++;
    }
  }
  return tparticles;
}

int64_t
generator_lattice_2d(
    const int lattice_type,
    const int domain_type,
    const point_t& bbox_min,
    const point_t& bbox_max,
    const double sph_sep,
    int64_t posid,
    bool count_only = true,
    double x[] = NULL,
    double y[] = NULL,
    double z[] = NULL)
{
   // Determine radius of the domain as a radius of inscribed sphere
   double radius = 0.0;
   if(domain_type == 1 || domain_type == 0) {
      radius = bbox_max[0] - bbox_min[0];
      if (radius > bbox_max[1] - bbox_min[1])
          radius = bbox_max[1] - bbox_min[1];
   }
   radius = 0.5*radius;

   // Central coordinates: in most cases this should be centered at (0,0)
   double x_c = (bbox_max[0]+bbox_min[0])/2.;
   double y_c = (bbox_max[1]+bbox_min[1])/2.;

   // Coordinate extents
   double xmin = bbox_min[0], xmax = bbox_max[0];
   double ymin = bbox_min[1], ymax = bbox_max[1];

   // Lattice spacing in x- and y-directions
   double dx = sph_sep;
   double dy = sph_sep*sqrt(3.)/2.;

   // True number of particles to be determined and returned
   int64_t tparticles = 0;

   if (lattice_type==0) { // rectangular lattice
     for(double y_p=ymin; y_p<ymax; y_p+=dx)
     for(double x_p=xmin; x_p<xmax; x_p+=dx)
     if(in_domain_2d(x_p,y_p,
                     x_c,y_c,
                     bbox_min,bbox_max,
                     radius,domain_type)){
       if(!count_only){
         x[posid] = x_p;
         y[posid] = y_p;
         z[posid] = 0.0;
       }
       tparticles++;
       posid++;
     } // if in domain
   }
   else { // triangular lattice
     for(double y_p=ymin,    yo=0; y_p<ymax; y_p+=dy,yo=1-yo)
     for(double x_p=xmin +yo*dx/2; x_p<xmax; x_p+=dx)
     if(in_domain_2d(x_p,y_p,
                     x_c,y_c,
                     bbox_min,bbox_max,
                     radius,domain_type)){
       if(!count_only){
         x[posid] = x_p;
         y[posid] = y_p;
         z[posid] = 0.0;
       }
       tparticles++;
       posid++;
     } // if in domain
   } // lattice
   return tparticles;
}


int64_t
generator_lattice_3d(
    const int lattice_type,
    const int domain_type,
    const point_t& bbox_min,
    const point_t& bbox_max,
    const double sph_sep,
    int64_t posid,
    bool count_only = true,
    double x[] = NULL,
    double y[] = NULL,
    double z[] = NULL)
{
   // Determine radius of the domain as a radius of inscribed sphere
   double radius = 0.0;
   if(domain_type==1 || domain_type==0) {
      radius = bbox_max[0] - bbox_min[0];
      for(size_t j=1;j<gdimension;j++)
        if (radius > bbox_max[j] - bbox_min[j])
          radius = bbox_max[j] - bbox_min[j];
   }
   radius = 0.5*radius;

   // Central coordinates: in most cases this should be centered at (0,0,0)
   double x_c = (bbox_max[0]+bbox_min[0])/2.;
   double y_c = (bbox_max[1]+bbox_min[1])/2.;
   double z_c = (bbox_max[2]+bbox_min[2])/2.;

   // True number of particles to be determined and returned
   int64_t tparticles = 0;

   // Coordinate extents
   double xmin = bbox_min[0], xmax = bbox_max[0];
   double ymin = bbox_min[1], ymax = bbox_max[1];
   double zmin = bbox_min[2], zmax = bbox_max[2];

   // Lattice spacing
   double dx = sph_sep;
   double dy = sph_sep*sqrt(3.)/2.;
   double dz = sph_sep*sqrt(2./3.);

   // The loop for lattice_type==0 (rectangular)
   if(lattice_type==0){
     for(double z_p=xmin; z_p<zmax; z_p+=dx)
     for(double y_p=ymin; y_p<ymax; y_p+=dx)
     for(double x_p=zmin; x_p<xmax; x_p+=dx)
     if(in_domain_3d(x_p,y_p,z_p,
                     x_c,y_c,z_c,
                     bbox_min,bbox_max,
                     radius,domain_type)) {
       tparticles++;
       if(!count_only){
         x[posid] = x_p;
         y[posid] = y_p;
         z[posid] = z_p;
       }
       posid++;
     } // if in domain
   }
   else if(lattice_type==1){//hcp lattice in 3D
     for(double z_p=zmin,         zo=0; z_p<zmax; z_p+=dz, zo=1-zo)
     for(double y_p=ymin-zo*dy/3, yo=0; y_p<ymax; y_p+=dy, yo=1-yo)
     for(double x_p=xmin+(yo-zo)*dx/2.; x_p<xmax; x_p+=dx)
     if(in_domain_3d(x_p,y_p,z_p,
                     x_c,y_c,z_c,
                     bbox_min,bbox_max,
                     radius,domain_type)){
       if(!count_only){
         x[posid] = x_p;
         y[posid] = y_p;
         z[posid] = z_p;
       }
       tparticles++;
       posid++;
     } // if in domain
   }
   else if(lattice_type==2) {//fcc lattice in 3D
     for(double z_p=zmin,         zl=0; z_p<zmax; z_p+=dz, zl=(zl+1)*(zl<3))
     for(double y_p=ymin-zl*dy/3, yo=0; y_p<ymax; y_p+=dy, yo=1-yo)
     for(double x_p=xmin+(yo-zl)*dx/2.; x_p<xmax; x_p+=dx)
     if(in_domain_3d(x_p,y_p,z_p,
                     x_c,y_c,z_c,
                     bbox_min,bbox_max,
                     radius,domain_type)) {
       if(!count_only){
         x[posid] = x_p;
         y[posid] = y_p;
         z[posid] = z_p;
       }
       tparticles++;
       posid++;
     } // if in domain
   } // lattice_type

   return tparticles;
}

// wrappers (because function pointers don't accept default parameters)
int64_t generate_lattice_1d(const int lattice_type, const int domain_type,
    const point_t& bbox_min, const point_t& bbox_max, const double sph_sep,
    int64_t posid, double * x, double * y, double * z) {
  return generator_lattice_1d(lattice_type,domain_type,
         bbox_min,bbox_max,sph_sep,posid, false, x,y,z);
}
int64_t count_lattice_1d(const int lattice_type, const int domain_type,
    const point_t& bbox_min, const point_t& bbox_max, const double sph_sep,
    int64_t posid) {
  return generator_lattice_1d(lattice_type,domain_type,
         bbox_min,bbox_max,sph_sep,posid);
}

int64_t generate_lattice_2d(const int lattice_type, const int domain_type,
    const point_t& bbox_min, const point_t& bbox_max, const double sph_sep,
    int64_t posid, double * x, double * y, double * z) {
  return generator_lattice_2d(lattice_type,domain_type,
         bbox_min,bbox_max,sph_sep,posid, false, x,y,z);
}
int64_t count_lattice_2d(const int lattice_type, const int domain_type,
    const point_t& bbox_min, const point_t& bbox_max, const double sph_sep,
    int64_t posid) {
  return generator_lattice_2d(lattice_type,domain_type,
         bbox_min,bbox_max,sph_sep,posid);
}

int64_t generate_lattice_3d(const int lattice_type, const int domain_type,
    const point_t& bbox_min, const point_t& bbox_max, const double sph_sep,
    int64_t posid, double * x, double * y, double * z) {
  return generator_lattice_3d(lattice_type,domain_type,
         bbox_min,bbox_max,sph_sep,posid, false, x,y,z);
}
int64_t count_lattice_3d(const int lattice_type, const int domain_type,
    const point_t& bbox_min, const point_t& bbox_max, const double sph_sep,
    int64_t posid) {
  return generator_lattice_3d(lattice_type,domain_type,
         bbox_min,bbox_max,sph_sep,posid);
}


// pointer types
typedef int64_t (*lattice_generate_function_t)(const int, const int,
    const point_t&, const point_t&, const double, int64_t,
    double*, double*, double*);
typedef int64_t (*particle_count_function_t)(const int, const int,
    const point_t&, const point_t&, const double, int64_t);
lattice_generate_function_t generate;
particle_count_function_t   count;

/**
 * @brief  Installs the 'generate' and 'count' function pointers
 */
void select() {
  switch(gdimension) {
  case 1:
    generate = generate_lattice_1d;
    count = count_lattice_1d;
    break;
  case 2:
    generate = generate_lattice_2d;
    count = count_lattice_2d;
    break;
  case 3:
    generate = generate_lattice_3d;
    count = count_lattice_3d;
    break;
  default:
    std::cerr << "you should not be here" << std::endl;
  }
}

} // namespace particle_lattice

