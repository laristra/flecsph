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
 *  domain_type  - int value, 0:box, 1:sphere
 *  bbox_min     - bounding box: bottom-left corner (or a bounding cube if domain)
 *  bbox_max     - bounding box: top-right corner   (is spherical                )
 *  sph_sep      - particle separation within the lattice
 *  posid        - enter the particle ID number from which the function is called
 *                 it is important for when the function is assigning position
 *                 coordinates to the particles
 *  count_only   - boolean for determing particle number (returned for allocating
 *                 proper arrays) or writing positions to those arrays
 *  x,y,z[]      - the arrays to be filled for the positions of each particle
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
 *             bbox_min    - minimum position for the total domain
 *             bbox_max    - maximum position for the total domain
 *             domain_type - int value, 0:cube, 1:sphere
 */
bool in_domain_2d(const double x, const double y, 
    const point_t& bbox_min, const point_t& bbox_max, const int domain_type) {
  // within_domain checks to see if the position is within the bounding box
  const double xmin = bbox_min[0], xmax = bbox_max[0],
               ymin = bbox_min[1], ymax = bbox_max[1];
  bool within_domain = x>=xmin && x<=xmax && y>=ymin && y<=ymax;

  // extra check for a circular domain
  if(domain_type==1) {
    const double r = 0.5*std::min(xmax-xmin, ymax-ymin);
    const double x0 = 0.5*(xmin + xmax);
    const double y0 = 0.5*(ymin + ymax);
    within_domain &= ((x-x0)*(x-x0)+(y-y0)*(y-y0) <r*r);
  }
  return within_domain;
}

bool in_domain_3d(const double x, const double y, const double z,
    const point_t& bbox_min, const point_t& bbox_max, const int domain_type) {
  // within_domain checks to see if the position is within the total domain
  const double xmin = bbox_min[0], xmax = bbox_max[0],
               ymin = bbox_min[1], ymax = bbox_max[1],
               zmin = bbox_min[2], zmax = bbox_max[2];
  bool within_domain = 
       x>=xmin && x<=xmax && y>=ymin && y<=ymax && z>=zmin && z<=zmax;

  // extra check for a spherical domain
  if(domain_type==1) {
    const double r = 0.5*std::min(std::min(xmax-xmin,ymax-ymin),(zmax-zmin));
    const double x0 = 0.5*(xmin + xmax);
    const double y0 = 0.5*(ymin + ymax);
    const double z0 = 0.5*(zmin + zmax);
    within_domain &= ((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0) < r*r);
  }
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
  // Central coordinates: in most cases this should be centered at 0
  double x_c = (bbox_max[0] + bbox_min[0])/2.;

  // True number of particles to be determined and returned
  int64_t tparticles = 0;

  // regular lattice in 1D
  double xmin = bbox_min[0], xmax = bbox_max[0];
  for(double x_p = xmin; x_p < xmax; x_p += sph_sep){
    tparticles++;
    if(!count_only){
      x[posid] = x_p;
      y[posid] = 0.0;
      z[posid] = 0.0;
    }
    posid++;
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
     if(in_domain_2d(x_p,y_p, bbox_min,bbox_max, domain_type)){
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
     if(in_domain_2d(x_p,y_p, bbox_min,bbox_max, domain_type)) {
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
     if(in_domain_3d(x_p,y_p,z_p, bbox_min,bbox_max, domain_type)) {
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
     if(in_domain_3d(x_p,y_p,z_p, bbox_min,bbox_max, domain_type)) {
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
     if(in_domain_3d(x_p,y_p,z_p, bbox_min,bbox_max, domain_type)) {
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

