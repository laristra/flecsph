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
 * @file density_profiles.h
 * @author Oleg Korobkin
 * @date January 2019
 * @brief Interface to select various density profiles
 * 
 * Notes: 
 *  - all spherical density profiles have support with radius R = 1;
 *  - all profiles are normalized to total mass M = 1.
 */

#ifndef DENSITY_PROFILES_H
#define DENSITY_PROFILES_H

#include <stdlib.h>
#include "user.h"
#include "tree.h"
#include <math.h>
#include <boost/algorithm/string.hpp>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))
namespace density_profiles {

  // spherical density profile function
  typedef double  (*radial_function_t)(const double);
  static radial_function_t spherical_density_profile = NULL;
  static radial_function_t spherical_mass_profile = NULL;
  static radial_function_t spherical_drho_dr = NULL;

  // constants for the mesa density
  static double mesa_rho0;
  static double mesa_q;        // ratio of the slope width to the radius

  // tabulated density profiles
  static std::vector<double> rad_grid;
  static std::vector<double> rho_grid;
  static std::vector<double> mass_grid;
  static std::vector<double> drhodr_grid;

  /**
   * @brief  constant uniform density in a domain of radius R = 1,
   *         normalized such that the total mass M = 1
   * @param  r     - spherical radius
   */
  double rho_constant_density(const double r) {
    double rho = 0.0;
    if constexpr (gdimension == 1)
      rho = 0.5;
    
    if constexpr (gdimension == 2)
      rho = 1.0/M_PI;

    if constexpr (gdimension == 3)
      rho = 0.75/M_PI;
    return rho;
  }

  double mass_constant_density(const double r) {
    double mass = 0.0;
    if constexpr (gdimension == 1)
      mass = r;
    
    if constexpr (gdimension == 2)
      mass = SQ(r);

    if constexpr (gdimension == 3)
      mass = CU(r);

    return mass;
  }

  double drhodr_constant_density(const double r) {
    return 0.;
  }

  /**
   * @brief  parabolic density
   * @param  r     - spherical radius
   */
  double rho_parabolic_density(const double r) {
    double rho = 0.0;
    if constexpr (gdimension == 1)
      rho = 0.75*(1. - SQ(r));
    
    if constexpr (gdimension == 2)
      rho = 2./M_PI*(1. - SQ(r));
    
    if constexpr (gdimension == 3)
      rho = 15./(8.*M_PI)*(1. - SQ(r));

    return rho;
  }

  double mass_parabolic_density(const double r) {
    double mass = 0.0;
    if constexpr (gdimension == 1)
      mass = 0.5*r*(3.0 - SQ(r));
    
    if constexpr (gdimension == 2)
      mass = SQ(r)*(2. - SQ(r));
    
    if constexpr (gdimension == 3)
      mass = 0.5*CU(r)*(5. - 3.*SQ(r));

    return mass;
  }

  double drhodr_parabolic_density(const double r) {
    double drhodr = 0.0;
    if constexpr (gdimension == 1)
      drhodr = -1.5*r;
    
    if constexpr (gdimension == 2)
      drhodr = -4./M_PI*r;
    
    if constexpr (gdimension == 3)
      drhodr = -15./(4.*M_PI)*r;

    return drhodr;
  }

  /**
   * @brief  spherical "mesa" density: flat top and steep slopes
   *
   *           / rho0                      if r < r0;
   *           |
   * rho(r) = <  rho0 (1 - (r-r0)^2/dr^2)  if r0 < r < 1;
   *           |
   *           \ 0                         if r > 1.
   *
   * @param  r     - spherical radius
   */
  double mesa_mass_helper(const double r) { 
    const double dr = mesa_q, r0 = 1. - mesa_q;
    double mm = 0.0;

    if constexpr (gdimension == 2)
      mm = SQ(r)*(1. - (.5*SQ(r) + SQ(r0))/SQ(dr))
         + (.5*CU(r0) + 4.*CU(r))*r0/(3.*SQ(dr));

    if constexpr (gdimension == 3)
      mm = CU(r0)/3. + (CU(r)-CU(r0))/3.*(1. - SQ(r0)/SQ(dr))
                     + r0*(SQ(r)-SQ(r0))*(SQ(r)+SQ(r0))/(2.*SQ(dr))
                     - (SQ(r)*CU(r)-SQ(r0)*CU(r0))/(5*SQ(dr)); 
    return mm;
  }

  double rho_mesa_density(const double r) {
    double rho = 0.0;
    const double r0 = 1. - mesa_q;
    if (r < 1.-mesa_q) 
      rho = mesa_rho0;
    else if (r < .9999)
      rho = mesa_rho0*(1. - SQ(r-r0)/SQ(mesa_q));
    return rho;
  }

  double mass_mesa_density(const double r) {
    const double dr = mesa_q, r0 = 1. - mesa_q;
    double m = 0.0;
    if constexpr (gdimension == 1) {
      if (r < 1.-mesa_q) 
        m = 2.*mesa_rho0*r;
      else if (r - 1. < 1e-12)
        m = 2.*mesa_rho0*(r - CU(r-r0)/(3.*SQ(dr)));
    }
    if constexpr (gdimension == 2) {
      if (r < 1.-mesa_q) 
        m = M_PI * mesa_rho0 * SQ(r);
      else if (r - 1. < 1e-12)
        m = M_PI*mesa_rho0*mesa_mass_helper(r);
    }
    if constexpr (gdimension == 3) {
      if (r < 1.-mesa_q) 
        m = 4.*M_PI/3. * mesa_rho0 * CU(r);
      else if (r - 1. < 1e-12)
        m = 4.*M_PI*mesa_rho0*mesa_mass_helper(r);
    }
    return m;
  }

  double drhodr_mesa_density(const double r) {
    double drhodr = 0.0;
    const double r0 = 1. - mesa_q;
    if (r > 1.-mesa_q and r < .9999) 
      drhodr = -2.*mesa_rho0*(r-r0)/SQ(mesa_q);
    return drhodr;
  }

  /**
   * @brief  read the density input file
   * @param  ifname - 4-column ASCII file: 1:r 2:rho 3:m 4:drho/dr
   *                  possibly with a header with lines starting with '#'
   */
  void read_input_density_file(const char * ifname) {
    using namespace std;
    ifstream infile;
    string line;
    int ln, Nr;
    std::string::size_type sz1, sz2;
    double rad, rho, mass, drhodr;

    // attempt to open the file
    infile.open (ifname);
    if (!infile) {
      cerr << "ERROR: Unable to open density profile '"<<ifname<<"'" <<endl;
      exit(1);
    }

    // count the number of lines, skipping the header
    Nr = 0;
    while (std::getline(infile,line)) {
      if (line.find("#") != string::npos or
          line.find_first_not_of(' ') == string::npos)
        continue;
      ++Nr;
    }
    cout << "Read density profile "<<ifname<<" with " << Nr 
         << " data points." << endl;

    // allocate arrays
    rad_grid.resize(Nr);
    rho_grid.resize(Nr);
    mass_grid.resize(Nr);
    drhodr_grid.resize(Nr);

    // read in the values
    infile.clear();
    infile.seekg(0, ios::beg);
    for(ln=0; std::getline(infile,line); ) {
      if (line.find("#") != string::npos or
          line.find_first_not_of(' ') == string::npos)
        continue;
      rad_grid[ln] = std::stod (line, &sz1);
      rho_grid[ln] = std::stof (line.substr(sz1), &sz2);
      sz1 += sz2;
      mass_grid[ln] = std::stof (line.substr(sz1), &sz2);
      sz1 += sz2;
      drhodr_grid[ln] = std::stof (line.substr(sz1));

      ln++;
    }
  }

  /**
   * @brief   get index i such that xp[i] < x < xp[i+1] (binary search)
   * @param   x     - the value to localize;
   * @param   xp    - array of increasing values where to localize x.
   * @return  index i of an interval containing point x;
   *          '-1'   if x < x[0];
   *          'N-1'  if x > x[N-1] (N is the vector size).
   */
  int get_interval_index ( const double x, const std::vector<double>& xp) {
    const int N = xp.size();
    if (x*(1 + 1e-15) < xp[0])   return -1;
    if (x*(1 - 1e-15) > xp[N-1]) return N-1;
    int i, i1 = 0, i2 = N-1;
    while (i2 - i1 > 1) {
      i = (i1 + i2)/2;
      if (x < xp[i]) 
        i2 = i;
      else
        i1 = i;
    }
    return i1;
  }


  /**
   * @brief   cubic interpolation
   * @param   x        - the x-coordinate location where to interpolate
   * @param   xp       - the x-coordinates of data points: increasing
   * @param   yp       - the y-coordinates of data points
   * @return  if x is inside the range of xp, returns interpolated value;
   *          if x is outside the range of xp: returns zero
   */
  double cubic_interp( const double x, const std::vector<double>& xp,  
                                       const std::vector<double>& yp) {
    int i2 = get_interval_index(x, xp);
    const int N = xp.size();
    if (i2 == 0)    i2 = 1;
    if (i2 > N - 3) i2 = N - 3;
    int
      i1 = i2 - 1,
      i3 = i2 + 1,
      i4 = i2 + 2;

    double 
      xx1 = x - xp[i1],
      xx2 = x - xp[i2],
      xx3 = x - xp[i3],
      xx4 = x - xp[i4];

    return yp[i1]*xx2/(xx2-xx1)*xx3/(xx3-xx1)*xx4/(xx4-xx1)
         + yp[i2]*xx1/(xx1-xx2)*xx3/(xx3-xx2)*xx4/(xx4-xx2)
         + yp[i3]*xx1/(xx1-xx3)*xx2/(xx2-xx3)*xx4/(xx4-xx3)
         + yp[i4]*xx1/(xx1-xx4)*xx2/(xx2-xx4)*xx3/(xx3-xx4);
  }


  /**
   * @brief  density profile from file, specified by
   *         the parameter density_profile_input
   * @param  r     - spherical radius
   */
  double rho_from_input_file(const double r) {
    return cubic_interp(r, rad_grid, rho_grid);
  }

  double mass_from_input_file(const double r) {
    return cubic_interp(r, rad_grid, mass_grid);
  }

  double drhodr_from_input_file(const double r) {
    return cubic_interp(r, rad_grid, drhodr_grid);
  }

  /**
   * @brief      Density profile selector
   */
  void select() {
    using namespace param;
    if (boost::iequals(density_profile,"constant")) {
      spherical_density_profile = rho_constant_density;
      spherical_mass_profile = mass_constant_density;
      spherical_drho_dr = drhodr_constant_density;
    }
    else if (boost::iequals(density_profile,"parabolic")) {
      spherical_density_profile = rho_parabolic_density;
      spherical_mass_profile = mass_parabolic_density;
      spherical_drho_dr = drhodr_parabolic_density;
    }
    else if (boost::iequals(density_profile,"mesa")) {
      spherical_density_profile = rho_mesa_density;
      spherical_mass_profile = mass_mesa_density;
      spherical_drho_dr = drhodr_mesa_density;
      mesa_q = mesa_rim_width;
      if constexpr (gdimension == 1) 
        mesa_rho0 = .5/(1. - mesa_q/3.);

      if constexpr (gdimension == 2)
        mesa_rho0 = 1./(M_PI*mesa_mass_helper(1.));

      if constexpr (gdimension == 3)
        mesa_rho0 = 1./(4.*M_PI*mesa_mass_helper(1.));
    }
    else if (boost::iequals(density_profile,"from file")) {
      // read rho input file
      read_input_density_file(input_density_file);
      spherical_density_profile = rho_from_input_file;
      spherical_mass_profile = mass_from_input_file;
      spherical_drho_dr = drhodr_from_input_file;
    }
    else {
      clog(error) << "ERROR: wrong parameter in density_profiles";
      exit(2);
    }

  } // select()

} // namespace density_profiles

#undef SQ
#undef CU
#endif // DENSITY_PROFILES_H
