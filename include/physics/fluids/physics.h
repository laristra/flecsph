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

#ifndef _physics_fluids_h_
#define _physics_fluids_h_

#include <cmath>
#include <vector>

#include "tree.h"
#include "io.h"

#include "kernel.h"

namespace physics{

  double dt = 0.000001;
  double K = 3.5;
  double gravity_cste = 9.80665; 
  point_t gravity_force;
  double totaltime = 0.0;
  double rest_density = 1000.;
  double alpha = 0.;
  double beta = 0.;
  double local_gamma = 7.;
  double cs0 = 10.;
  double eta_sq = 0.01;
  double background_pressure = 0.;
  bool do_verlet_cor = false; 
  double coeff_dt = 0.3;
  int verlet_cstep = 10.;
  int verlet_current = 0.;
  double epsilon = 0.1;
  double maxtime = 10.;
  double outputtime = 0.02;
  double MAC = 0.;
  

  auto kernel = kernel::quintic_wendland_2D;
  auto kernel_gradient = kernel::quintic_wendland_2D_gradient;

  
  void init_physics(
    const char * filename)
  {
    verlet_cstep = io::input_parameter_int(filename,"verlet_cstep");
    eta_sq = io::input_parameter_double(filename,"eta_sq");
    cs0 = io::input_parameter_double(filename,"rest_sound_speed");
    K = io::input_parameter_double(filename,"gas_stiffness");
    gravity_cste = io::input_parameter_double(filename,"gravity_cste");
    dt = io::input_parameter_double(filename,"timestep");
    alpha = io::input_parameter_double(filename,"alpha");
    beta = io::input_parameter_double(filename,"beta");
    local_gamma = io::input_parameter_double(filename,"gamma");
    epsilon = io::input_parameter_double(filename,"epsilon");
    maxtime = io::input_parameter_double(filename,"maxtime");
    outputtime = io::input_parameter_double(filename,"outputtime");
    rest_density = io::input_parameter_double(filename,"rest_density");

    //maxtime = 0.5;
    // Switch kernel for 3D case
    if(gdimension == 3){
        kernel = kernel::quintic_wendland_3D;
        kernel_gradient = kernel::quintic_wendland_3D_gradient;
    }

    /*printf("\nInput Data:\n"
      "verlet_cstep=%d\n"
      "eta_sq=%g\n"
      "cs0=%g\n"
      "K=%g\n"
      "gravity_cste=%g\n"
      "dt=%g\n"
      "alpha=%g\n"
      "beta=%g\n"
      "gamma=%g\n"
      "epsilon=%g\n"
      "maxtime=%g\n"
      "outputtime=%g\n"
      "rest_density=%g\n\n",verlet_cstep,eta_sq,cs0,K,gravity_cste,dt,alpha,
      beta,local_gamma,epsilon,maxtime,outputtime,rest_density);
    */

    if(gdimension == 2 ){
      gravity_force[0] = 0.;
      gravity_force[1] = -gravity_cste;
    }else{
      gravity_force[0] = 0.;
      gravity_force[1] = -gravity_cste;
      gravity_force[2] = 0.;
    }
  }

  double mu_ij(
    double vijrij, 
    double r, 
    double h)
  {
    return h*vijrij/(r*r+eta_sq);
  }


  double viscosity(
    double vijrij, 
    double h, 
    double r,
    double rhoab, 
    double csab)
  {
    if(vijrij < 0){
      double muij = mu_ij(vijrij,r,h);
      return (-alpha*csab*muij+beta*muij*muij)/(rhoab);
    }
    return 0.; 
  }

  void 
  init_density_velocity(
    body_holder* srch)
  {
    body* source = srch->getBody();
    source->setVelocityTmp(source->getVelocity());
    source->setDensityTmp(source->getDensity());
    assert(source->getDensityTmp() > 0. && 
      source->getDensityTmp() == source->getDensity());
  }

  void 
  compute_accel(
      body_holder* srch,
      std::vector<body_holder*> nbsh)
  {
    body* source = srch->getBody();
    // For verlet 
    source->setVelocityNM1(source->getVelocityTmp());
    source->setVelocityTmp(source->getVelocity());

    assert(source->getDensityTmp() > 0.);
    assert(source->getDensity() > 0.);
    source->setDensityNM1(source->getDensityTmp());
    source->setDensityTmp(source->getDensity());
    
    point_t accel = {};
    double density_dt = 0.;
    point_t velCor = {};
    double max_visc = -1000.;
    for(auto nbh: nbsh){
      // All the neighbors are in the right distance already 
      body* nb = nbh->getBody();
      double r = flecsi::distance(source->getPosition(),nb->getPosition());
      if( 1.0e-18 < r && r <= 2*source->getSmoothinglength()){       
        
        point_t vij = source->getVelocity() - nb->getVelocity();
        point_t rij = source->getPosition() - nb->getPosition();
        double vijrij = flecsi::dot(
          flecsi::point_to_vector(vij),
          flecsi::point_to_vector(rij));

        // Compute the kernels
        double wab = kernel(r,source->getSmoothinglength());
        point_t grad_wab = kernel_gradient(rij,source->getSmoothinglength());
        
        double density_bar = (source->getDensity()+nb->getDensity())/2.;
        assert(!std::isnan(density_bar));
        double cs_bar = (source->getSoundspeed()+nb->getSoundspeed())/2.;
        assert(!std::isnan(cs_bar));

        // Pressure contribution
        double pressure_part = source->getPressure()/pow(source->getDensity(),2)+
          nb->getPressure()/pow(nb->getDensity(),2);
        double viscosity_part = viscosity(
          vijrij,source->getSmoothinglength(),r,density_bar,cs_bar);

        assert(!std::isnan(pressure_part));
        assert(!std::isnan(viscosity_part));

        accel = accel - nb->getMass()*(pressure_part+viscosity_part)*grad_wab;
        //printf("accel = %g %g with m=%g pv=%g with (di=%g dj=%g)"
        //  " vp=%g grad=%g %g\n",
        //  accel[0],accel[1],nb->getMass(),pressure_part,source->getDensity(),
        //  nb->getDensity(),viscosity_part,
        //  grad_wab[0],grad_wab[1]);
        //fflush(stdout);
        assert(!std::isnan(accel[0])&&!std::isnan(accel[1]));
        
        // Compute the velocity correction
        velCor = velCor - wab * vij * nb->getMass() / density_bar; 

        // Save the biggest local visc value 
        max_visc = std::max(max_visc,vijrij/(r*r+eta_sq));
        
        // Compute the denisty variation 
        density_dt = density_dt + nb->getMass()*
          flecsi::dot(
            flecsi::point_to_vector(vij),
            flecsi::point_to_vector(grad_wab));
      }
    }
    assert(!std::isnan(accel[0])&&!std::isnan(accel[1]));
    source->setAcceleration(accel+gravity_force);
    source->setVelocityCor(source->getVelocity()+epsilon*velCor);
    source->setDensityDt(density_dt);
    source->setMaxVisc(max_visc);
  } // compute_accel
 
  void 
  update_dt(
    std::vector<body_holder*> bodies)
  {
    // Consider the verlet 
    verlet_current++;
    do_verlet_cor = false;
    if(verlet_current % verlet_cstep == 0){
      do_verlet_cor = true;
      verlet_current = 0.;
    }

    double h = bodies[0]->getBody()->getSmoothinglength();
    // Search for max soundspeed
    double cs_max = -1000.;
    double visc_max = -1000.;
    double dt1 = sqrt(h/sqrt(
      flecsi::dot(
        flecsi::point_to_vector(gravity_force),
        flecsi::point_to_vector(gravity_force))));
    double dt2 = 1;
    // Compute the new dt 
    for(auto b: bodies){
      // Compute the Fa_max
      dt1 = std::min(dt1,
        sqrt(h/sqrt(flecsi::dot(
          flecsi::point_to_vector(
            b->getBody()->getAcceleration()+gravity_force),
          flecsi::point_to_vector(
            b->getBody()->getAcceleration()+gravity_force)))));
      cs_max = std::max(b->getBody()->getSoundspeed(),cs_max);
      visc_max = std::max(visc_max,b->getBody()->getMaxVisc());
    }

    dt2 = dt1;
    if(cs_max > 0.){
      dt2 = h/(cs_max+h*visc_max);
    }

    dt = std::min(dt1,dt2)*coeff_dt;
  }

  void 
  EOS(
      body_holder* srch)
  { 
    body* source = srch->getBody();
    double pressure = K*(pow(source->getDensity()/rest_density,local_gamma)-1.)
      + background_pressure;
    double soundspeed = cs0*pow(source->getDensity()/rest_density,3);
    assert(!std::isnan(pressure) && !std::isnan(soundspeed));
    source->setPressure(pressure);
    source->setSoundspeed(soundspeed);
  } // EOS 

  void 
  verlet_integration_EOS(
      body_holder* srch)
  {
    body* source = srch->getBody();

    point_t position = source->getPosition();
    point_t velocity = {};
    double density = 0.;
    assert(source->getDensityNM1()>0.);

    if(!do_verlet_cor){
      if(!source->is_wall()){
        position += source->getVelocityCor()*dt+
          source->getAcceleration()*dt*dt*0.5;
        velocity = source->getVelocityNM1()+source->getAcceleration()*2*dt;
      }
      density = source->getDensityNM1() + source->getDensityDt()*2.*dt;
    }else{
      if(!source->is_wall()){
        position += source->getVelocityCor()*dt+
          source->getAcceleration()*dt*dt*0.5;
        velocity = source->getVelocity()+source->getAcceleration()*dt;
      }
      density = source->getDensity() + source->getDensityDt()*dt;
    }
    for(size_t i = 0; i < gdimension; ++i){
      assert(!std::isnan(velocity[i]&&!std::isnan(velocity[i])));
      assert(!std::isnan(position[i]&&!std::isnan(position[i])));
    }
    assert(!std::isnan(density));
    assert(density > 0.);

    source->setVelocity(velocity);
    source->setPosition(position);
    source->setDensity(density);

    // Compute EOS
    EOS(srch); 
  } // verlet_integration_EOS
   
}; // fluid

#endif // _physics_physics_h_
