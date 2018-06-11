/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@ @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@///// /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@      /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@ /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////  /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@      /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@      /@@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
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
  double gamma = 2.0; // 1.4; 
  double K = 1;
  double epsilon = 1;
  double g_strength = 1; 
  double damp = 1;
  double totaltime = 0.0;
  double MAC = 0.;
  double eta = 0.01;
  double A = 0.6366197723675814;
  double angular_moment = 1.36049047255;
  //double rotation_speed = 10.;
  double angular_speed = 0.48;
  double QZZ = 0.;
  double Omega2 = 0.;
  double t_relax = 1.12;

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
      double kernelresult = (1./2.)*(
          kernel(dist,source->getSmoothinglength())+
          kernel(dist,nb->getSmoothinglength()));
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
    double pressure = source->getAdiabatic()*
      pow(source->getDensity(),gamma);
    source->setPressure(pressure);
  } // compute_pressure


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
    double soundspeed = pow(gamma*
        source->getPressure()/source->getDensity(),1./2.);
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
      body* nb,
      double gamma)
  {  
    double result = 0.0;
    double h_ij = (1./2.)*
      (source->getSmoothinglength()+nb->getSmoothinglength()); 
    space_vector_t vecVelocity = flecsi::point_to_vector(
        source->getVelocity() - nb->getVelocity());
    space_vector_t vecPosition = flecsi::point_to_vector(
        source->getPosition() - nb->getPosition());
    double dotproduct = flecsi::dot(vecVelocity,vecPosition);
    if(dotproduct >= 0.0)
      return result;
    // Should add norm to space_vector
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    result = dotproduct / (h_ij*(((dist*dist)/(h_ij*h_ij))+eta*eta));
    mpi_assert(result < 0.0);
    return result; 
  } // mu

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
    point_t acceleration = source->getAcceleration();

    // Compute the adiabatic factor here 
    double dadt = 0;
    point_t hydro = {};
    for(auto nbh : ngbsh){ 
      body* nb = nbh->getBody();

      if(nb->getPosition() == source->getPosition()){
        continue;
      }

      // Artificial viscosity
      double density_ij = (1./2.)*(source->getDensity()+nb->getDensity());
      double soundspeed_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
      double mu_ij = mu(source,nb,epsilon);
      double viscosity = (-alpha*mu_ij*soundspeed_ij+beta*mu_ij*mu_ij)
        /density_ij;
      mpi_assert(viscosity>=0.0);


      // Hydro force
      double pressureDensity = 
          source->getPressure()/(source->getDensity()*source->getDensity())
          + nb->getPressure()/(nb->getDensity()*nb->getDensity());

      point_t vecPosition = source->getPosition()-nb->getPosition();
      point_t sourcekernelgradient = gradKernel(
          vecPosition,source->getSmoothinglength());
      point_t nbkernelgradient = gradKernel(
          vecPosition,nb->getSmoothinglength());
      point_t resultkernelgradient = (1./2.)*
        (sourcekernelgradient+nbkernelgradient);

      hydro += nb->getMass()*(pressureDensity+viscosity)
        *resultkernelgradient;

      // Compute the adiabatic factor evolution 
      dadt += nb->getMass() * viscosity * 
        flecsi::dot(
          flecsi::point_to_vector(source->getVelocity()-nb->getVelocity()),
          flecsi::point_to_vector(resultkernelgradient)
        );
    }
    hydro = -1.0*hydro;
    acceleration += hydro;
    source->setAcceleration(acceleration);
    dadt *= (gamma - 1)/(2*pow(source->getDensity(),gamma-1));
    source->setDadt(dadt);

    // Integrate the Ai factor 

  } // compute_hydro_acceleration


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

  void 
  apply_rotation(
    body_holder* srch)
  {
    body* source = srch->getBody();
    double theta = angular_speed * dt;
    mpi_assert(theta > 0);
    point_t tmp = source->getPosition();
    point_t new_position = tmp;
    new_position[0] = tmp[0]*cos(theta)-tmp[1]*sin(theta);
    new_position[1] = tmp[0]*sin(theta)+tmp[1]*cos(theta);
    source->setPosition(new_position);
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
      double local_mu = fabs(mu(source,nb,gamma));
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
  compute_QZZ(
    std::vector<body_holder*>& bodies)
  {
    QZZ = 0.;
    Omega2 = 0.;
    for(auto bh: bodies){
      body * src = bh->getBody();
      QZZ += src->getMass() * 
        (src->getPosition()[0]*src->getPosition()[0]+
          src->getPosition()[1]*src->getPosition()[1]);
    }
  }

  void compute_rotation(
    std::vector<body_holder*>& bodies)
  {
    Omega2 = pow(angular_moment/QZZ,2);

    for(auto bh: bodies){
      body * src = bh->getBody();
      point_t frot = src->getPosition();
      frot[0] *= Omega2;
      frot[1] *= Omega2;
      if(gdimension == 3){
        frot[2] = 0.;
      }
      src->setAcceleration(src->getAcceleration() + frot);
      // Correct the velocity 
      src->setAcceleration(src->getAcceleration() - 
        src->getVelocity()/t_relax);
    }
  }

  // Internal energy from Adiabatic factor 
  // u = A / (g-1) * rho ^ {g - 1}
  void 
  compute_internal_energy(
    body_holder* srch)
  {
    body* source = srch->getBody();
    double u = source->getAdiabatic()/(gamma-1)*
      pow(source->getDensity(),gamma-1);
    source->setInternalenergy(u);
  }

  // Adiabatic ratio from internal energy 
  void 
  compute_adiabatic_from_u(
    body_holder* srch)
  {
    body* source = srch->getBody();
    double A = (gamma-1)*source->getInternalenergy()/
      pow(source->getDensity(),gamma-1);
    source->setAdiabatic(A);
  }

}; // physics

#endif // _physics_physics_h_
