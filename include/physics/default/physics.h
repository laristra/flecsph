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

#ifndef _physics_physics_h_
#define _physics_physics_h_

#include <vector>

#include "kernel.h"
#include "tree.h"

namespace physics{

  bool fixed_timestep = true;
  bool do_boundaries = false;
  bool stop_boundaries = false; 
  bool reflect_boundaries = false;
  point_t max_boundary = {};
  point_t min_boundary = {};
  double dt = 0.0;
  double alpha = 2; 
  double beta = 2; // 1; 
  double gamma = 1.4;//2.0; // 1.4; 
  double K = 1;
  double epsilon = 1;
  double g_strength = 1; 
  double damp = 1;
  double totaltime = 0.0;
  double MAC = 0.;
  double eta = 0.01;

  // Default configuration for kernel
  int kernel_choice = 0;
  auto kernel = kernel::cubic_spline_kernel;
  auto gradKernel = kernel::cubic_spline_gradKernel;

  // Compute density based on neighbors
  // The formula used here is:
  // rho_i = \sum_j m_j*cubic\_spline\_kernel\_3D_ij
  //static
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

  // Computre pressure on a body 
  // This function does not need the neighbors
  // Formula is:
  // P_i = u_i*rho_i^{\Gamma}
  //static
  void 
  compute_pressure(
      body_holder* srch)
  { 
    body* source = srch->getBody();
    double pressure = (gamma-1.0)*
      source->getDensity()*source->getInternalenergy();
    source->setPressure(pressure);
  } // compute_pressure

  //For zero temperature white dwarf EOS
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

  // Compute the sound speed on a body 
  // This function does not need the neighbors
  // Formula is:
  // CS_i = (P_i/rho_i)^{1/2}
  //static
  void 
  compute_soundspeed(
      body_holder* srch)
  {
    body* source = srch->getBody();
    double soundspeed = sqrt(gamma*source->getPressure()/source->getDensity());
    source->setSoundspeed(soundspeed);
  } // computeSoundspeed

  void 
  compute_density_pressure_soundspeed(
    body_holder* srch, 
    std::vector<body_holder*>& nbsh)
  {
    compute_density(srch,nbsh);
    compute_pressure(srch);
    compute_soundspeed(srch); 
  }

  // mu used in artificial viscosity 
  // Formula is:
  // h_ij*(v_i-v_j).(p_i-p_j)/(dist*dist + h_ij+epsi)
  // with h_ij = (h_i+h_j)/2
  //static
  double 
  mu(
      body* source, 
      body* nb)
  {  
    double result = 0.0;
    double h_ij = .5*(source->getSmoothinglength()+nb->getSmoothinglength()); 
    space_vector_t vecVelocity = flecsi::point_to_vector(
        source->getVelocity() - nb->getVelocity());
    space_vector_t vecPosition = flecsi::point_to_vector(
        source->getPosition() - nb->getPosition());
    double dotproduct = flecsi::dot(vecVelocity,vecPosition);

    if(dotproduct >= 0.0)
      return result;
    // Should add norm to space_vector
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    result = h_ij * dotproduct / (dist*dist + epsilon);
    
    mpi_assert(result < 0.0);
    return result; 
  } // mu

  double 
  viscosity(
    body* source, 
    body* nb)
  {
    double rho_ij = (1./2.)*(source->getDensity()+nb->getDensity());
    double c_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
    double mu_ij = mu(source,nb);
    double res = (-alpha*c_ij*mu_ij+beta*mu_ij*mu_ij)/rho_ij;
    mpi_assert(res>=0.0);
    return res;
  }

  // Compute the hydro force with the artificial viscosity 
  // Formula is:
  // f_{hydro}_i = m_j * (\sum_j P_i/rho_i^2 + P_j/rho_j^2 + Pi_{ij})*
  // gradKernel_{ij}
  // with Pi_{ij} the artificial viscosity 
  //static
  void 
  compute_hydro_acceleration(
    body_holder* srch, 
    std::vector<body_holder*>& ngbsh)
  { 
    body* source = srch->getBody();

    // Add in the acceleration
    point_t acceleration = point_t{};//source->getAcceleration();

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
      point_t nbkernelgradient = gradKernel(
          vecPosition,nb->getSmoothinglength());
      point_t resultkernelgradient = (1./2.)*
        (sourcekernelgradient+nbkernelgradient);

      hydro += nb->getMass()*(pressureDensity+visc)
        *resultkernelgradient;

    }
    hydro = -1.0*hydro;
    acceleration += hydro;
    source->setAcceleration(acceleration);
    //std::cout<<*source<<std::endl;

  } // compute_hydro_acceleration

  // compute u 
  // Formula is:
  // u_i = P_i/rho_i^2 \sum_j m_j*v_ij*gradKernel_ij
  //static
  void
  compute_internalenergy(
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
      point_t nbkernelgradient = gradKernel(
          vecPosition,nb->getSmoothinglength());
      space_vector_t resultkernelgradient = flecsi::point_to_vector((1./2.)*
          (sourcekernelgradient+nbkernelgradient));

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
  } // compute_internalenergy


  void dudt_integration(
      body_holder* srch)
  {
    body* source = srch->getBody(); 
    source->setInternalenergy(
      source->getInternalenergy()+dt*source->getDudt());
  }

  // Apply boundaries 
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
      for(int i = 0; i < gdimension; ++i){
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
    //std::cout<<"dt1 = "<<dt1<<std::endl;
  
    // Second based on max mu 
    double max_mu_ij = -9999999;
    for(auto nbh: ngbhs){
      body* nb = nbh->getBody(); 
      double local_mu = fabs(mu(source,nb));
      max_mu_ij = std::max(max_mu_ij,local_mu);
    }
    double dt2 = source->getSmoothinglength()/
        (source->getSoundspeed()+
         1.2*alpha*source->getSoundspeed()+
         1.2*beta*max_mu_ij);
    dt2 *= 0.1;

    //std::cout<<"dt2 = "<<dt1<<std::endl;

    double min = std::min(dt1,dt2);
  #pragma omp critical 
    physics::dt = std::min(dt,min);
    //return std::min(physics::dt,dt1); 
  }

  // Calculate simple linear momentum for checking momentum conservation
  void 
  compute_lin_momentum(
      std::vector<body_holder*>& bodies, 
      point_t* total) 
  {
    *total = {0};
    for(auto nbh: bodies) {
      *total += nbh->getBody()->getLinMomentum(); 
    }  
  }

  void
  nsquare_gravitation(
      std::vector<body_holder*>& bodies){
    int64_t nelem = bodies.size();
    #pragma omp parallel for 
    for(int64_t i = 0; i<nelem; ++i){
      body * source = bodies[i]->getBody();
      point_t grav = point_t{};
      for(auto& n: bodies){
        body* nb = n->getBody();
        double dist = flecsi::distance(source->getPosition(),nb->getPosition());
        point_t pos = source->getPosition() - nb->getPosition();
        if(dist > 0.){
          grav += -nb->getMass() / (dist*dist*dist) * pos; 
        }
      }
      source->setAcceleration(source->getAcceleration() + grav);
    }
  }

}; // physics

#endif // _physics_physics_h_
