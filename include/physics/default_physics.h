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

// Basic elements for the physics
namespace physics{
  double dt = 0.0;
  double totaltime = 0.0;
  int64_t iteration = 0;
}

#include "params.h"
#include "utils.h"
#include "kernels.h"
#include "tree.h"
#include "eforce.h"

#include "eos.h"
#include "integration.h"
#include "boundary.h"
#include "viscosity.h"

namespace physics{
  using namespace param;

  /**
   * @brief      Compute the density
   * Based on Fryer/05 eq(10)
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void
  compute_density(
      body* source,
      std::vector<body*>& nbsh)
  {
    double density = 0;
    const int n_nb = nbsh.size();
    mpi_assert(nbsh.size()>0);
    // Array of distances
    std::vector<double> distances(n_nb);
    std::vector<double> masses(n_nb);
    std::vector<double> nb_radius(n_nb);
    double radius = source->radius();
    point_t coordinates = source->coordinates();

    for(int i = 0 ; i < n_nb; ++i){
      masses[i] = nbsh[i]->mass();
      nb_radius[i] = nbsh[i]->radius();
      distances[i] = flecsi::distance(coordinates,nbsh[i]->coordinates());
    }
    for(int i = 0 ; i < n_nb; ++i){
      double kernel = 0;
      if constexpr (gdimension == 1){
        kernel = kernels::wendland_c6_1d(
          distances[i],.5*(radius+nb_radius[i]));
      }else{
        kernel = kernels::wendland_c6_23d(
          distances[i],.5*(radius+nb_radius[i]));
      }
      density += kernel*masses[i];
    } // for
    mpi_assert(density>0);
    source->setDensity(density);
  } // compute_density


  /**
   * @brief      Calculates total energy for every particle
   * @param      srch  The source's body holder
   */
  void set_total_energy (body_holder* srch) {
    body* source = srch->getBody();
    const point_t pos = source->coordinates(),
                  vel = source->getVelocity();
    const double eint = source->getInternalenergy(),
                 epot = external_force::potential(pos);
    double ekin = vel[0]*vel[0];
    for (unsigned short i=1; i<gdimension; ++i)
      ekin += vel[i]*vel[i];
    ekin *= .5;
    source->setTotalenergy(eint + epot + ekin);
  } // set_total_energy


  /**
   * @brief      Subtracts mechanical energy from total energy
   *             to recover internal energy
   * @param      srch  The source's body holder
   */
  void recover_internal_energy (body* source) {
    const point_t pos = source->coordinates(),
                  vel = source->getVelocity();
    const double etot = source->getTotalenergy(),
                 epot = external_force::potential(pos);
    double ekin = vel[0]*vel[0];
    for (unsigned short i=1; i<gdimension; ++i)
      ekin += vel[i]*vel[i];
    ekin *= .5;
    const double eint = etot - ekin - epot;
    if (eint < 0.0) {
      std::cerr << "ERROR: internal energy is negative!" << std::endl
                << "particle id: " << source->id()    << std::endl
                << "total energy: " << etot              << std::endl
                << "kinetic energy: " << ekin            << std::endl
                << "particle position: " << pos          << std::endl;
      mpi_assert(false);
    }
    source->setInternalenergy(eint);
  } // recover_internal_energy


  /**
   * @brief      Compute the density, EOS and spundspeed in the same function
   * reduce time to gather the neighbors
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void
  compute_density_pressure_soundspeed(
    body* source,
    std::vector<body*>& nbsh)
  {
    compute_density(source,nbsh);
    if (thermokinetic_formulation)
      recover_internal_energy(source);
    eos::compute_pressure(source);
    eos::compute_soundspeed(source);
  }


  /**
   * @brief      Calculates the hydro acceleration
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void
  compute_acceleration(
    body* source,
    std::vector<body*>& ngbsh)
  {
    using namespace param;
    // Reset the accelerastion
    // \TODO add a function to reset in main_driver
    point_t acceleration = {};
    point_t hydro = {};
    const double h_s = source->radius();
    const point_t coordinates = source->coordinates();
    const double density = source->getDensity();
    const double pressure = source->getPressure();
    source->setMumax(0.0);
    const int n_nb = ngbsh.size();

    double viscosities[n_nb];
    point_t sourcekernelgradient[n_nb];

    double densities[n_nb];
    double pressures[n_nb];
    point_t positions[n_nb];
    double radii[n_nb];
    double masses[n_nb];

    int index=-1;
    for(int i = 0 ; i < n_nb; ++i){
      densities[i] = ngbsh[i]->getDensity();
      pressures[i] = ngbsh[i]->getPressure();
      positions[i] = ngbsh[i]->coordinates();
      radii[i] = ngbsh[i]->radius();
      masses[i] = ngbsh[i]->mass();
      if(positions[i]==coordinates)
        index=i;
    }

    for(int i = 0 ; i < n_nb; ++i){ // not vectorized
      viscosities[i] = viscosity::viscosity(source,ngbsh[i]);
    }
    for(int i = 0 ; i < n_nb; ++i){ // Not vectorized due to Kernel
      // Kernel computation
      point_t vecPosition = coordinates - positions[i];
      if constexpr (gdimension == 1 ){
        sourcekernelgradient[i] = kernels::gradient_wendland_c6_1d(
          vecPosition,(h_s+radii[i])*.5);
      }else{
        sourcekernelgradient[i] = kernels::gradient_wendland_c6_23d(
          vecPosition,(h_s+radii[i])*.5);
      }
    } // for 
    //ignore itself
    sourcekernelgradient[index]={};
    viscosities[index]=0;
    for(int i = 0 ; i < n_nb; ++i){ // Vectorized
      double pressureDensity = pressure/(density*density)
          + pressures[i]/(densities[i]*densities[i]);
      hydro += masses[i]*(pressureDensity + viscosities[i])
        *sourcekernelgradient[i];
    }
    hydro = -1.0*hydro;
    acceleration += hydro;
    acceleration += external_force::acceleration(source);
    source->setAcceleration(acceleration);
  } // compute_hydro_acceleration


  /**
   * @brief      Adds drag force to acceleration
   * @param      srch  The source's body holder
   */
  void add_drag_acceleration( body_holder* srch) {
    using namespace param;
    body* source = srch->getBody();
    point_t       acc = source->getAcceleration();
    const point_t vel = source->getVelocity();
    acc += external_force::acceleration_drag(vel);
    source->setAcceleration(acc);
  } // add_drag_acceleration


  /**
   * @brief      Calculates the dudt, time derivative of internal energy.
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void compute_dudt(
      body* source,
      std::vector<body*>& ngbsh)
  {
    double dudt = 0;
    double dudt_pressure = 0.;
    double dudt_visc = 0.;
    const double h_s = source->radius();

    for(auto nb: ngbsh){
      if(nb->coordinates() == source->coordinates()){
        continue;
      }

      // Artificial viscosity
      double visc = viscosity::viscosity(source,nb);

      // Compute the gradKernel ij
      const double h_n = nb->radius();
      point_t vecPosition = source->coordinates()-nb->coordinates();
      point_t sourcekernelgradient = kernels::gradKernel(
          vecPosition,(h_s+h_n)*.5);
      space_vector_t resultkernelgradient =
          flecsi::point_to_vector(sourcekernelgradient);

      // Velocity vector
      space_vector_t vecVelocity = flecsi::point_to_vector(
          source->getVelocity() - nb->getVelocity());

      dudt_pressure += nb->mass()*
        flecsi::dot(vecVelocity,resultkernelgradient);
      dudt_visc += visc*nb->mass()*
        flecsi::dot(vecVelocity,resultkernelgradient);
    }

    double P_a = source->getPressure();
    double rho_a = source->getDensity();
    dudt = P_a/(rho_a*rho_a)*dudt_pressure + .5*dudt_visc;

    // Do not change internal energy in relaxation phase
    if(iteration < relaxation_steps){
       dudt = 0.0;
    }

    source->setDudt(dudt);
  } // compute_dudt


  /**
   * @brief      Calculates the dedt, time derivative of either
   *             thermokinetic (internal + kinetic) or total
   *             (internal + kinetic + potential) energy.
   * See e.g. Rosswog (2009) "Astrophysical SPH" eq. (34)
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void compute_dedt(
      body* source,
      std::vector<body*>& ngbsh)
  {

    double dedt = 0;

    const point_t pos_a = source->coordinates(),
                  vel_a = source->getVelocity();
    const double h_a = source->radius(),
                 P_a = source->getPressure(),
                 rho_a = source->getDensity();
    const double Prho2_a = P_a/(rho_a*rho_a);

    const int n_nb = ngbsh.size();
    point_t pos_b[n_nb];
    double h_b[n_nb];
    point_t vel_b[n_nb];
    double rho_b[n_nb];
    double P_b[n_nb];
    double m_b[n_nb];

    point_t Da_Wab[n_nb];
    double Pi_ab[n_nb];


    int index=-1;
    for(int i = 0 ; i < n_nb; ++i){
      rho_b[i] = ngbsh[i]->getDensity();
      P_b[i] = ngbsh[i]->getPressure();
      pos_b[i] = ngbsh[i]->coordinates();
      h_b[i] = ngbsh[i]->radius();
      vel_b[i] = ngbsh[i]->getVelocity();
      m_b[i] = ngbsh[i]->mass();
      if(pos_b[i]==pos_a)
        index=i;
    }

    for(int i = 0 ; i < n_nb; ++i){
      // Compute the \nabla_a W_ab
      Da_Wab[i] = kernels::gradKernel(pos_a - pos_b[i], (h_a+h_b[i])*.5);
      Pi_ab[i] = viscosity::viscosity(source,ngbsh[i]);
    }

    // To avoid the same particle
    m_b[index] = 0;
    Da_Wab[index] = {};
    Pi_ab[index] = 0;

    for(int i = 0 ; i < n_nb; ++i){ // Vectorized
      // va*DaWab and vb*DaWab
      point_t tmp_a = vel_a*Da_Wab[i];
      point_t tmp_b = vel_b[i]*Da_Wab[i];
      double va_dot_DaWab = tmp_a[0];
      double vb_dot_DaWab = tmp_b[0];
      for(unsigned int j = 1 ; j < gdimension ; ++j){
        va_dot_DaWab += tmp_a[j];
        vb_dot_DaWab += tmp_b[j];
      }
      const double Prho2_b = P_b[i]/(rho_b[i]*rho_b[i]);
      // add this neighbour's contribution
      dedt -= m_b[i]*( Prho2_a*vb_dot_DaWab + Prho2_b*va_dot_DaWab
                + .5*Pi_ab[i]*(va_dot_DaWab + vb_dot_DaWab));
    }
    source->setDedt(dedt);
  } // compute_dedt



  /**
   * @brief      Adds energy dissipation rate due to artificial
   *             particle relaxation drag force
   * @param      srch  The source's body holder
   */
  void add_drag_dedt( body_holder* srch) {
    using namespace param;
    body* source = srch->getBody();
    const point_t vel = source->getVelocity();
    const point_t acc = external_force::acceleration_drag(vel);
    double va = vel[0]*acc[0];
    for (short int i=1; i<gdimension; ++i)
      va += vel[i]*acc[i];
    double dedt = source->getDedt();
    source->setDedt(dedt + va);
  } // add_drag_dedt


  /**
   * @brief      Compute the timestep from acceleration and mu
   * From CES-Seminar 13/14 - Smoothed Particle Hydrodynamics
   *
   * @param      srch   The source's body holder
   */
  void compute_dt(body_holder* srch) {
    body* source = srch->getBody();
    const double tiny = 1e-24;
    const double mc   = 0.6; // constant in denominator for viscosity

    // particles separation around this particle
    const double dx = source->radius()
                    / (sph_eta*kernels::kernel_width);

    // timestep based on particle velocity
    const point_t vel = source->getVelocity();
    const double vn  = norm_point(vel);
    const double dt_v = dx/(vn + tiny);

    // timestep based on acceleration
    const double acc = norm_point(source->getAcceleration());
    const double dt_a = sqrt(dx/(acc + tiny));

    // timestep based on sound speed and viscosity
    const double max_mu_ab = source->getMumax();
    const double cs_a = source->getSoundspeed();
    const double dt_c = dx/ (tiny + cs_a*(1 + mc*sph_viscosity_alpha)
                                  + mc*sph_viscosity_beta*max_mu_ab);

    // minimum timestep
    double dtmin = timestep_cfl_factor * std::min(std::min(dt_v,dt_a), dt_c);

    // timestep based on positivity of internal energy
    if (thermokinetic_formulation) {
      const double eint = source->getInternalenergy();
      const point_t pos = source->coordinates();
      const double epot = external_force::potential(pos);
      double epot_next;
      int i;
      for(i=0; i<20; ++i) {
        epot_next = external_force::potential(pos + dtmin*vel);
        if(epot_next - epot < eint*0.5) break;
        dtmin *= 0.5;
      }
      assert (i<20);
    }

    source->setDt(dtmin);
  }


  /**
   * @brief      Reduce adaptive timestep and set its value
   *
   * @param      bodies   Set of bodies
   */
  void set_adaptive_timestep(
      std::vector<body_holder>& bodies
      )
  {
    double dtmin = 1e24; // some ludicrous number

    #pragma omp parallel for reduction(min:dtmin)
    for(size_t i = 0 ; i < bodies.size(); ++i){
      if(!bodies[i].is_local()) continue;
      dtmin = std::min(dtmin, bodies[i].getBody()->getDt());
    }

    mpi_utils::reduce_min(dtmin);

    if (dtmin < physics::dt)
      physics::dt = std::min(dtmin, physics::dt/2.0);

    if (dtmin > 2.0*physics::dt)
      physics::dt = physics::dt*2.0;
  }


  void
  compute_smoothinglength(
      std::vector<body_holder>& bodies)
  {
    if (gdimension == 1) {
      #pragma omp parallel for
      for(size_t i = 0 ; i < bodies.size(); ++i){
        if(!bodies[i].is_local()) continue;
        auto particle = bodies[i].getBody();
        double m_b   = particle->mass();
        double rho_b = particle->getDensity();
        particle->set_radius(
          m_b/rho_b * sph_eta*kernels::kernel_width);
      }
    }
    else if (gdimension == 2) {
      #pragma omp parallel for
      for(size_t i = 0 ; i < bodies.size(); ++i){
        if(!bodies[i].is_local()) continue;
        auto particle = bodies[i].getBody();
        double m_b   = particle->mass();
        double rho_b = particle->getDensity();
        particle->set_radius(
          sqrt(m_b/rho_b) * sph_eta*kernels::kernel_width);
      }
    }
    else {
      #pragma omp parallel for
      for(size_t i = 0 ; i < bodies.size(); ++i){
        if(!bodies[i].is_local()) continue;
        auto particle = bodies[i].getBody();
        double m_b   = particle->mass();
        double rho_b = particle->getDensity();
        particle->set_radius(
          cbrt(m_b/rho_b) * sph_eta*kernels::kernel_width);
      }
    } // if gdimension
  }


  /**
   * @brief update smoothing length for particles (Rosswog'09, eq.51)
   *
   * ha = eta/N \sum_b pow(m_b / rho_b,1/dimension)
   */
  void compute_average_smoothinglength(
      std::vector<body_holder>& bodies,
      int64_t nparticles) {
    compute_smoothinglength(bodies);
    // Compute the total
    double total = 0.;

    #pragma omp parallel for reduction(+:total)
    for(size_t i = 0 ; i < bodies.size(); ++i)
    {
      if(!bodies[i].is_local()) continue;
      total += bodies[i].getBody()->radius();
    }

    // Add up with all the processes
    MPI_Allreduce(MPI_IN_PLACE,&total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // Compute the new smoothing length
    double new_h = 1./(double)nparticles * total;
    #pragma omp parallel for
    for(size_t i = 0 ; i < bodies.size(); ++i){
      if(!bodies[i].is_local()) continue;
      bodies[i].getBody()->set_radius(new_h);
    }
  }
}; // physics

#endif // _default_physics_h_
