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
 * @file physics.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic physics implementation
 */

#ifndef _default_physics_h_
#define _default_physics_h_

#include <vector>

#include "params.h"
#include "utils.h"
#include "kernels.h"
#include "tree.h"
#include <boost/algorithm/string.hpp>

namespace physics{
  using namespace param;

  /**
   * Default values
   * They are usually modified in the main_driver 
   * to be application-specific. 
   * \TODO add a parameter file to be read in the main_driver
   */
  point_t max_boundary = {};
  point_t min_boundary = {};
  double dt = 0.0;
  double damp = 1;
  double totaltime = 0.0;
  double A = 1.0;
  double MAC = 1.;
  int64_t iteration = 0;

  /**
   * @brief      Kernel selector: types, global variables and the function
   *
   * @param      kstr     Kernel string descriptor
   *
   * @return     Pointer to the kernel
   */
  typedef double  (*kernel_function_t)(const double,    const double);
  typedef point_t (*kernel_gradient_t)(const point_t &, const double);
  kernel_function_t kernel;
  kernel_gradient_t gradKernel;

  void
  select_kernel(const std::string& kstr) {
    if (boost::iequals(kstr,"cubic spline")) {
      kernel = kernels::cubic_spline;
      gradKernel = kernels::gradient_cubic_spline;
    }
    else if (boost::iequals(kstr, "quintic spline")) {
      kernel = kernels::quintic_spline;
      gradKernel = kernels::gradient_quintic_spline;
    }
    else if (boost::iequals(kstr, "gaussian")) {
      kernel = kernels::gaussian;
      gradKernel = kernels::gradient_gaussian;
    }
    else if (boost::iequals(kstr, "wendland quintic")) {
      kernel = kernels::wendland_quintic;
      gradKernel = kernels::gradient_wendland_quintic;
    }
    else {
      clog_one(fatal) << "Bad kernel parameter" << std::endl;
    }
  }


  /**
   * @brief      Compute the density 
   * Based on Fryer/05 eq(10)
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void 
  compute_density(
      body_holder* srch, 
      std::vector<body_holder*>& nbsh)
  {
    body* source = srch->getBody();
    double density = 0;
    mpi_assert(nbsh.size()>0);
    for(auto nbh : nbsh){
      body* nb = nbh->getBody();
      double dist = flecsi::distance(source->getPosition(),nb->getPosition());
      mpi_assert(dist>=0.0);
      double kernelresult = kernel(dist,
            .5*(source->getSmoothinglength()+nb->getSmoothinglength()));
      density += kernelresult*nb->getMass();
    } // for
    mpi_assert(density>0);
    source->setDensity(density);
  } // compute_density

  /**
   * @brief      Compute the pressure
   * Ideal gas EOS
   * @param      srch  The source's body holder
   */
  void 
  compute_pressure(
      body_holder* srch)
  { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = (poly_gamma-1.0)*
      source->getDensity()*source->getInternalenergy();
    source->setPressure(pressure);
  } // compute_pressure

#ifdef ADIABATIC
  /**
   * @brief      Compute the pressure based on adiabatic index
   *
   * @param      srch  The source's body holder
   */
  void 
  compute_pressure_adiabatic(
      body_holder* srch)
  { 
    using namespace param;
    body* source = srch->getBody();
    double pressure = source->getAdiabatic()*
      pow(source->getDensity(),poly_gamma);
    source->setPressure(pressure);
  } // compute_pressure
#endif 

  /**
   * @brief      \TODO NEED COMMENTS
   *
   * @param      srch  The srch
   */
  void 
  compute_pressure_wd(
      body_holder* srch)
  { 
    body* source = srch->getBody();
    double A_dwd = 6.00288e22;
    double B_dwd = 9.81011e5;

    double x_dwd = pow((source->getDensity())/B_dwd,1.0/3.0);
    double pressure = A_dwd*(x_dwd*(2.0*x_dwd*x_dwd-3.0)*
 		      pow(x_dwd*x_dwd+1.0,1.0/2.0)+3.0*asinh(x_dwd));
    source->setPressure(pressure);
  } // compute_pressure_wd



  /**
   * @brief      Compute the sound speed
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics 
   * 
   * @param      srch  The source's body holder
   */
  void 
  compute_soundspeed(
      body_holder* srch)
  {
    using namespace param;
    body* source = srch->getBody();
    double soundspeed = sqrt(poly_gamma*source->getPressure()
                                       /source->getDensity());
    source->setSoundspeed(soundspeed);
  } // computeSoundspeed

  /**
   * @brief      Compute the density, EOS and spundspeed in the same function 
   * reduce time to gather the neighbors
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void 
  compute_density_pressure_soundspeed(
    body_holder* srch, 
    std::vector<body_holder*>& nbsh)
  {
    compute_density(srch,nbsh);
    compute_pressure(srch);
    compute_soundspeed(srch); 
  }

#ifdef ADIABATIC
  void 
  compute_density_pressure_adiabatic_soundspeed(
    body_holder* srch, 
    std::vector<body_holder*>& nbsh)
  {
    compute_density(srch,nbsh);
    compute_pressure_adiabatic(srch);
    compute_soundspeed(srch); 
  }
#endif

  /**
   * @brief      mu_ij for the artificial viscosity 
   * From Rosswog'09 (arXiv:0903.5075) - 
   * Astrophysical Smoothed Particle Hydrodynamics, eq.(60) 
   *
   * @param      srch     The source particle
   * @param      nbsh     The neighbor particle
   *
   * @return     Contribution for mu_ij of this neighbor
   *
   * @uses       epsilon  global parameter
   */
  double 
  mu(
      body* source, 
      body* nb)
  {  
    using namespace param;
    double result = 0.0;
    double h_ij = .5*(source->getSmoothinglength()+nb->getSmoothinglength()); 
    space_vector_t vecVelocity = flecsi::point_to_vector(
        source->getVelocity() - nb->getVelocity());
    space_vector_t vecPosition = flecsi::point_to_vector(
        source->getPosition() - nb->getPosition());
    double dotproduct = flecsi::dot(vecVelocity,vecPosition);

    if(dotproduct >= 0.0)
      return result;
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    result = h_ij*dotproduct / (dist*dist + sph_viscosity_epsilon*h_ij*h_ij);
    
    mpi_assert(result < 0.0);
    return result; 
  } // mu

  /**
   * @brief      Artificial viscosity 
   * From Rosswog'09 (arXiv:0903.5075) - 
   * Astrophysical Smoothed Particle Hydrodynamics, eq.(59) 
   *
   * @param      srch  The source particle
   * @param      nbsh  The neighbor particle
   *
   * @return     The artificial viscosity contribution 
   */
  double 
  viscosity(
    body* source, 
    body* nb)
  {
    using namespace param;
    double rho_ij = (1./2.)*(source->getDensity()+nb->getDensity());
    double c_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
    double mu_ij = mu(source,nb);
    double res = ( -sph_viscosity_alpha*c_ij*mu_ij
                  + sph_viscosity_beta*mu_ij*mu_ij)/rho_ij;
    mpi_assert(res>=0.0);
    return res;
  }

  /**
   * @brief      Calculates the hydro acceleration
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics 
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void 
  compute_hydro_acceleration(
    body_holder* srch, 
    std::vector<body_holder*>& ngbsh)
  { 
    body* source = srch->getBody();

    // Reset the accelerastion 
    // \TODO add a function to reset in main_driver
    point_t acceleration = point_t{};

    point_t hydro = {};
    for(auto nbh : ngbsh){ 
      body* nb = nbh->getBody();

      if(nb->getPosition() == source->getPosition()){
        continue;
      }

      // Compute viscosity
      double visc = viscosity(source,nb);
      
      // Hydro force
      point_t vecPosition = source->getPosition()-nb->getPosition();
      double pressureDensity = source->getPressure()/(source->getDensity()*
          source->getDensity())
          + nb->getPressure()/(nb->getDensity()*nb->getDensity());

      // Kernel computation
      point_t sourcekernelgradient = gradKernel(
          vecPosition,source->getSmoothinglength());
      point_t resultkernelgradient = sourcekernelgradient;

      hydro += nb->getMass()*(pressureDensity+visc)
        *resultkernelgradient;

    }
    hydro = -1.0*hydro;
    acceleration += hydro;
    source->setAcceleration(acceleration);
  } // compute_hydro_acceleration

  /**
   * @brief      Calculates the dudt, variation of internal energy.
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics 
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void
  compute_dudt(
    body_holder* srch,
    std::vector<body_holder*>& ngbsh)
  {
    body* source = srch->getBody();

    double dudt = 0;
    double dudt_pressure = 0.;
    double dudt_visc = 0.;

    for(auto nbh: ngbsh){
      body* nb = nbh->getBody();

      if(nb->getPosition() == source->getPosition()){
        continue;
      }

      // Artificial viscosity
      double visc = viscosity(source,nb);
    
      // Compute the gradKernel ij      
      point_t vecPosition = source->getPosition()-nb->getPosition();
      point_t sourcekernelgradient = gradKernel(
          vecPosition,source->getSmoothinglength());
      space_vector_t resultkernelgradient = 
          flecsi::point_to_vector(sourcekernelgradient);

      // Velocity vector 
      space_vector_t vecVelocity = flecsi::point_to_vector(
          source->getVelocity()-nb->getVelocity());

      dudt_pressure += nb->getMass()*
        flecsi::dot(vecVelocity,resultkernelgradient);
      dudt_visc += visc*nb->getMass()*
        flecsi::dot(vecVelocity,resultkernelgradient);
    }
    
    dudt = source->getPressure()/(source->getDensity()*source->getDensity())*
    dudt_pressure+.5*dudt_visc;

    source->setDudt(dudt);
  } // compute_dudt



#ifdef ADIABATIC
  /**
   * @brief      Compute the adiabatic index for the particles 
   *
   * @param      srch   The source's body holder
   * @param      ngbsh  The neighbors' body holders
   */
  void 
  compute_dadt(
    body_holder* srch, 
    std::vector<body_holder*>& ngbsh)
  { 
    using namespace param;
    body* source = srch->getBody();
    
    // Compute the adiabatic factor here 
    double dadt = 0;
    
    for(auto nbh : ngbsh){ 
      body* nb = nbh->getBody();

      if(nb->getPosition() == source->getPosition()){
        continue;
      }

      // Artificial viscosity
      double density_ij = (1./2.)*(source->getDensity()+nb->getDensity());
      double soundspeed_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
      double mu_ij = mu(source,nb);
      double viscosity = ( -sph_viscosity_alpha*mu_ij*soundspeed_ij
                          + sph_viscosity_beta*mu_ij*mu_ij)/density_ij;
      mpi_assert(viscosity>=0.0);

      point_t vecPosition = source->getPosition()-nb->getPosition();
      point_t sourcekernelgradient = gradKernel(
          vecPosition,source->getSmoothinglength());
      point_t resultkernelgradient = sourcekernelgradient;

      
      // Compute the adiabatic factor evolution 
      dadt += nb->getMass() * viscosity * 
        flecsi::dot(
          flecsi::point_to_vector(source->getVelocity()-nb->getVelocity()),
          flecsi::point_to_vector(resultkernelgradient)
        );
    }
    
    dadt *= (poly_gamma - 1) / 
            (2*pow(source->getDensity(),poly_gamma-1));
    source->setDadt(dadt);
 

  } // compute_hydro_acceleration
#endif 

  /**
   * @brief      Integrate the internal energy variation, update internal energy
   *
   * @param      srch  The source's body holder
   */
  void dudt_integration(
      body_holder* srch)
  {
    body* source = srch->getBody(); 
    source->setInternalenergy(
      source->getInternalenergy()+dt*source->getDudt());
  }

#ifdef ADIABATIC 
    /**
   * @brief      Integrate the internal energy variation, update internal energy
   *
   * @param      srch  The source's body holder
   */
  void dadt_integration(
      body_holder* srch)
  {
    body* source = srch->getBody(); 
    source->setAdiabatic(
      source->getAdiabatic()+dt*source->getDadt());
  }
#endif 

  /**
   * @brief      Apply boundaries if they are set
   *
   * @param      srch  The source's body holder
   *
   * @return     True if the particle have been considered outside the 
   * boundaries
   */
  bool
  compute_boundaries(
      body_holder* srch)
  {
    body* source = srch->getBody();
    point_t velocity = source->getVelocity();
    point_t position = source->getPosition();
    point_t velocityHalf = source->getVelocityhalf();

    bool considered = false;

    if(stop_boundaries){
      bool stop = false; 
      for(size_t i = 0; i < gdimension; ++i){
        if(position[i] < min_boundary[i] ||
          position[i] > max_boundary[i]){
          stop = true; 
        }
      }
      if(stop){
        velocity = point_t{};
        velocityHalf = point_t{};
        considered = true;
      
      }
    }else if(reflect_boundaries){
      for(size_t dim=0;dim < gdimension ; ++dim){
        if(position[dim] < min_boundary[dim] || 
            position[dim] > max_boundary[dim]){
          double barrier = max_boundary[dim];
          if(position[dim] < min_boundary[dim]){
            barrier = min_boundary[dim];
          }

          // Here just invert the velocity vector and velocityHalf 
          double tbounce = (position[dim]-barrier)/velocity[dim];
          position -= velocity*(1-damp)*tbounce;

          position[dim] = 2*barrier-position[dim];
          velocity[dim] = -velocity[dim];
          velocityHalf[dim] = -velocityHalf[dim];

          velocity *= damp;
          velocityHalf *= damp;
          considered = true;
        }
      }
    }
    source->setPosition(position);
    source->setVelocity(velocity);
    source->setVelocityhalf(velocityHalf);
    return considered;
  }

#if 0 
  // \TODO VERSION USED IN THE BNS, CHECK VALIDITY REGARDING THE OTHER ONE 
  void 
  leapfrog_integration(
      body_holder* srch)
  {
    body* source = srch->getBody();
    

    point_t velocity = source->getVelocityhalf()+
      source->getAcceleration() * dt / 2.;
    point_t velocityHalf = velocity+
      source->getAcceleration() * dt / 2.;
    point_t position = source->getPosition()+velocityHalf*dt;
    // integrate dadt 
    double adiabatic_factor = source->getAdiabatic() + source->getDadt()* dt;

    source->setVelocity(velocity);
    source->setVelocityhalf(velocityHalf);
    source->setPosition(position);
    source->setAdiabatic(adiabatic_factor);
    
    mpi_assert(!std::isnan(position[0])); 
  }
#endif

  /**
   * @brief      Leapfrog integration, first step 
   *
   * @param      srch  The source's body holder
   */
  void 
  leapfrog_integration_first_step(
      body_holder* srch)
  {
    body* source = srch->getBody();

    // If wall, reset velocity and dont move 
    if(source->is_wall()){
      source->setVelocity(point_t{});
      source->setVelocityhalf(point_t{}); 
      return; 
    }

    point_t velocityHalf = source->getVelocity() + 
        dt/2.*source->getAcceleration();
    point_t position = source->getPosition()+velocityHalf*dt;
    point_t velocity = 1./2.*(source->getVelocityhalf()+velocityHalf);

    if(do_boundaries){
      if(physics::compute_boundaries(srch)){
        return;
      }
    }

    source->setVelocityhalf(velocityHalf);
    source->setVelocity(velocity);
    source->setPosition(position);

    mpi_assert(!std::isnan(position[0])); 
  }

  /**
   * @brief      Leapfrog integration
   *
   * @param      srch  The source's body holder
   */
  void 
  leapfrog_integration(
      body_holder* srch)
  {
    body* source = srch->getBody();
    
    // If wall, reset velocity and dont move 
    if(source->is_wall()){
      source->setVelocity(point_t{});
      source->setVelocityhalf(point_t{}); 
      return; 
    }
    
    point_t velocityHalf = source->getVelocityhalf() + 
        dt*source->getAcceleration();
    point_t position = source->getPosition()+velocityHalf*dt;
    point_t velocity = 1./2.*(source->getVelocityhalf()+velocityHalf);

    if(do_boundaries){
      if(physics::compute_boundaries(srch)){
        return;
      }
    }

    source->setVelocityhalf(velocityHalf);
    source->setVelocity(velocity);
    source->setPosition(position);
    
    mpi_assert(!std::isnan(position[0])); 
  }

  /**
   * @brief      Compute the timestep from acceleration and mu 
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics 
   *
   * @param      srch   The source's body holder
   * @param      ngbhs  The neighbors' body holders
   */
  void 
  compute_dt(
      body_holder* srch,
      std::vector<body_holder*>& ngbhs)
  {
    body* source = srch->getBody();
    
    // First compute dt based on acceleration norm 
    double accelNorm = 0.0;
    for(size_t i=0;i<gdimension;++i){
      accelNorm += source->getAcceleration()[i]*source->getAcceleration()[i];
    }
    accelNorm = sqrt(accelNorm);
    double dt1 = pow(source->getSmoothinglength()/accelNorm,1.0/2.0);
  
    // Second based on max mu 
    double max_mu_ij = -9999999;
    for(auto nbh: ngbhs){
      body* nb = nbh->getBody(); 
      double local_mu = fabs(mu(source,nb));
      max_mu_ij = std::max(max_mu_ij,local_mu);
    }
    double dt2 = source->getSmoothinglength()/
        (source->getSoundspeed()+
         1.2*sph_viscosity_alpha*source->getSoundspeed()+
         1.2*sph_viscosity_beta*max_mu_ij);
    dt2 *= 0.1;

    double min = std::min(dt1,dt2);
    // Critical OMP to avoid outside synchronizations
  #pragma omp critical 
    physics::dt = std::min(dt,min);
  }




}; // physics

#endif // _default_physics_h_
