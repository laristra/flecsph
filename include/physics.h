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

#ifndef _physics_h_
#define _physics_h_

#include <vector>

#include "kernel.h"
#include "tree.h"

namespace physics{

  bool fixed_timestep = true;
  double dt = 0.0;
  double alpha = 2; 
  double beta = 1; 
  double gamma = 2; 
  double epsilon = 1; 

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
    assert(nbsh.size()>0);
    for(auto nbh : nbsh){
      body* nb = nbh->getBody();
      double dist = flecsi::distance(source->getPosition(),nb->getPosition());
      assert(dist>=0.0);
      double kernelresult = (1./2.)*(
          kernel(dist,source->getSmoothinglength())+
          kernel(dist,nb->getSmoothinglength()));
      assert(kernelresult>=0.0);
      density += kernelresult*nb->getMass();
    } // for
    assert(density>0);
    source->setDensity(density);
  } // c ompute_density

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
    source->setPressure(source->getInternalenergy()*
      pow(source->getDensity(),gamma));
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
    source->setSoundspeed(pow(gamma*
        source->getPressure()/source->getDensity(),1./2.));
  } // computeSoundspeed

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
    result = h_ij * dotproduct / (dist*dist + epsilon*h_ij);
    assert(result < 0.0);
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

    point_t hydro = {};
    for(auto nbh : ngbsh){ 
      body* nb = nbh->getBody();

      // Artificial viscosity
      double density_ij = (1./2.)*(source->getDensity()+nb->getDensity());
      double soundspeed_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
      double mu_ij = mu(source,nb,epsilon);
      double viscosity = (-alpha*mu_ij*soundspeed_ij+beta*mu_ij*mu_ij)
        /density_ij;
      assert(viscosity>=0.0);

      // Hydro force
      point_t vecPosition = source->getPosition()-nb->getPosition();
      double pressureDensity = source->getPressure()/(source->getDensity()*
          source->getDensity())
          + nb->getPressure()/(nb->getDensity()*nb->getDensity());
      point_t sourcekernelgradient = gradKernel(
          vecPosition,source->getSmoothinglength());
      point_t nbkernelgradient = gradKernel(
          vecPosition,nb->getSmoothinglength());
      point_t resultkernelgradient = (1./2.)*
        (sourcekernelgradient+nbkernelgradient);

      hydro += nb->getMass()*(pressureDensity+viscosity)
        *resultkernelgradient;
    }
    hydro = -1.0*hydro;
    acceleration += hydro;
    source->setAcceleration(acceleration);
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
    double internalenergy = 0.0;

    for(auto nbh: ngbsh){
      body* nb = nbh->getBody();
    
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

      // Add artificial viscosity
      // \TODO compute it one for acceleration and internalenergy 
      double density_ij = (1./2.)*(source->getDensity()+nb->getDensity());
      double soundspeed_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
      double mu_ij = mu(source,nb,epsilon);
      double viscosity = (-alpha*mu_ij*soundspeed_ij+beta*mu_ij*mu_ij)
        /density_ij;
      assert(viscosity>=0.0);

      internalenergy += nb->getMass()*
        (source->getPressure()/(source->getDensity()*source->getDensity())
        + 1./2.*viscosity)*
        flecsi::dot(vecVelocity,resultkernelgradient);
    }

    //internalenergy *= source->getPressure()/
    //  (source->getDensity()*source->getDensity());
    source->setInternalenergy(source->getInternalenergy()+internalenergy);
  } // compute_internalenergy

  // Compute acceleration for the FGrav and FHydro 
  // If we consider the relaxation step, we need FRoche or FRot 
  //static 
  void 
  compute_acceleration(
      body_holder* srch)
  {
    body* source = srch->getBody();
    point_t acceleration = {0.,0.,0.};
    acceleration += source->getHydroForce() / source->getMass();
    acceleration += source->getGravForce();
    //acceleration /= source->getMass();
  
    // Use it in the new velocity 
    point_t velocity = {0.,0.,0.};
    velocity = source->getVelocityhalf() + dt/2.0*acceleration;

    source->setAcceleration(acceleration);
    source->setVelocity(velocity);
  }

  // Leapfrog integration to move the particles 
  // velHalf = vel+dt/2*acc
  // pos = pos + dt*velHalf
  //static
  void 
  leapfrog_integration(
    body_holder* srch)
  {
    body* source = srch->getBody();
    point_t position;
    point_t velocityhalf; 
    if(source->getPosition()[0] > 0.1 && source->getPosition()[0]<1){
      velocityhalf = source->getVelocity()+dt/2.0*source->getAcceleration();
      position = source->getPosition()+ dt*velocityhalf ;
      source->setPosition(position);
    }
    //source->setVelocityhalf(velocityhalf);  
  } // leapfrog_integration  

}; // physics

#endif // _physics_h_
