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
#include "user.h"
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
   * @brief      Computes the density in "vanilla sph" formulation
   *             [Rosswog'09, eq.(13)]:
   *             
   *             $\rho_a =\sum_b {m_b W_ab(r_ab, (h_a + h_b)/2)}$
   *
   * @param      particle  The particle body 
   * @param      nbs       Vector of neighbor particles
   */
  void
  compute_density(
      body* particle,
      std::vector<body*>& nbs)
  {
    using namespace kernels;
    const double h_a = particle->radius();
    const point_t pos_a = particle->coordinates();
    const int n_nb = nbs.size();
    mpi_assert(n_nb>0);

    double r_a_[n_nb], m_[n_nb], h_[n_nb];
    for(int b = 0 ; b < n_nb; ++b){
      const body * const nb = nbs[b];
      m_[b]  = nb->mass();
      h_[b]  = nb->radius();
      point_t pos_b = nb->coordinates();
      r_a_[b] = flecsi::distance(pos_a, pos_b);
    }

    double rho_a = 0.0;
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double Wab =  kernel<kernel_type,gdimension>(r_a_[b],.5*(h_a+h_[b]));
      rho_a += m_[b]*Wab;
    } // for
    mpi_assert(rho_a>0);
    particle->setDensity(rho_a);
  } // compute_density


  /**
   * @brief      Calculates total energy for every particle
   * @param      srch  The source's body holder
   */
  void set_total_energy (body_holder* srch) {
    body* particle = srch->getBody();
    const point_t pos = particle->coordinates(),
                  vel = particle->getVelocity();
    const double eint = particle->getInternalenergy(),
                 epot = external_force::potential(pos);
    //const space_vector_t & svel = *reinterpret_cast<const space_vector_t *> (&vel);
    //double ekin = flecsi::dot(vel,vel)/2.0;
    double ekin = vel[0]*vel[0];
    for (unsigned short i=1; i<gdimension; ++i)
      ekin += vel[i]*vel[i];
    ekin *= .5;
    particle->setTotalenergy(eint + epot + ekin);
  } // set_total_energy


  /**
   * @brief      Subtracts mechanical energy from total energy
   *             to recover internal energy
   * @param      srch  The source's body holder
   */
  void recover_internal_energy (body* particle) {
    const point_t pos = particle->coordinates(),
                  vel = particle->getVelocity();
    const double etot = particle->getTotalenergy(),
                 epot = external_force::potential(pos);
    double ekin = vel[0]*vel[0];
    for (unsigned short i=1; i<gdimension; ++i)
      ekin += vel[i]*vel[i];
    ekin *= .5;
    const double eint = etot - ekin - epot;
    if (eint < 0.0) {
      std::cerr << "ERROR: internal energy is negative!" << std::endl
                << "particle id: " << particle->id()     << std::endl
                << "total energy: " << etot              << std::endl
                << "kinetic energy: " << ekin            << std::endl
                << "particle position: " << pos          << std::endl;
      mpi_assert(false);
    }
    particle->setInternalenergy(eint);
  } // recover_internal_energy


  /**
   * @brief      Compute the density, EOS and spundspeed in the same function
   * reduce time to gather the neighbors
   *
   * @param      particle  The particle body 
   * @param      nbs       Vector of neighbor particles
   */
  void
  compute_density_pressure_soundspeed(
    body* particle,
    std::vector<body*>& nbs)
  {
    compute_density(particle,nbs);
    if (thermokinetic_formulation)
      recover_internal_energy(particle);
    eos::compute_pressure(particle);
    eos::compute_soundspeed(particle);
  }


  /**
   * @brief      Calculates the hydro acceleration ("vanilla ice")
   *             [Rosswog'09, eqs.(29,55)]:
   *
   *     (dv_a)                (  P_a       P_b            )
   *     (----)   = -sum_b m_b ( -----  +  -----   + Pi_ab ) D_i Wab  + ext_i
   *     ( dt )_i              (rho_a^2   rho_b^2          )
   * 
   * @param      particle  The particle body 
   * @param      nbs       Vector of neighbor particles
   */
  void
  compute_acceleration(
    body* particle,
    std::vector<body*>& nbs)
  {
    using namespace param;
    using namespace viscosity;
    using namespace kernels;

    // Reset the accelerastion
    // \TODO add a function to reset in main_driver

    // this particle (index 'a')
    const double h_a = particle->radius(),
               rho_a = particle->getDensity(),
                 P_a = particle->getPressure(),
                 c_a = particle->getSoundspeed();
    const point_t pos_a = particle->coordinates(),
                  vel_a = particle->getVelocity();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double rho_[n_nb],P_[n_nb],h_[n_nb],m_[n_nb],c_[n_nb],Pi_a_[n_nb];
    point_t pos_[n_nb], vel_[n_nb], DiWa_[n_nb];

    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      rho_[b] = nb->getDensity();
      P_[b]   = nb->getPressure();
      pos_[b] = nb->coordinates();
      vel_[b] = nb->getVelocity();
      c_[b]   = nb->getSoundspeed();
      h_[b]   = nb->radius();
      m_[b]   = nb->mass() * (pos_[b]!=pos_a); // if same particle, m_b->0
    }

    // precompute viscosity and kernel gradients
    particle->setMumax(0.0);  // needed for adaptive timestep calculation
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double mu_ab = mu(h_a,h_[b],vel_a,vel_[b],pos_a,pos_[b]);
      Pi_a_[b] = artificial_viscosity(rho_a,rho_[b],c_a,c_[b],mu_ab);
      point_t pos_ab = pos_a - pos_[b];
      DiWa_[b] = 
          kernel_gradient<kernel_type,gdimension>(pos_ab, .5*(h_a + h_[b]));
    }

    // compute the final answer
    const double Prho2_a = P_a/(rho_a*rho_a);
    point_t acc_a = 0.0;
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      const double Prho2_b = P_[b]/(rho_[b]*rho_[b]);
      acc_a += -m_[b]*(Prho2_a + Prho2_b + Pi_a_[b])*DiWa_[b];
    }
    acc_a += external_force::acceleration(particle);
    particle->setAcceleration(acc_a);
  } // compute_hydro_acceleration


  /**
   * @brief      Adds drag force to acceleration
   * @param      srch  The source's body holder
   */
  void add_drag_acceleration( body_holder* srch) {
    using namespace param;
    body* particle = srch->getBody();
    point_t       acc = particle->getAcceleration();
    const point_t vel = particle->getVelocity();
    acc += external_force::acceleration_drag(vel);
    particle->setAcceleration(acc);
  } // add_drag_acceleration


  /**
   * @brief      Calculates the dudt, time derivative of internal energy.
   *             [Rosswog'09, eqs.(29,55)]:
   *
   *             du_a             (  P_a      1       )
   *             ---- = sum_b m_b ( -----  +  - Pi_ab ) (D_i Wab . v_ab)
   *              dt              (rho_a^2    2       )
   * 
   * @param      particle  The particle body 
   * @param      nbs       Vector of neighbor particles
   */
  void compute_dudt(
      body* particle,
      std::vector<body*>& nbs)
  {
    // Do not change internal energy in relaxation phase
    if(iteration < relaxation_steps){
       particle->setDudt(0.0);
       return;
    }

    using namespace viscosity;
    using namespace kernels;

    // this particle (index 'a')
    const double h_a = particle->radius(),
               rho_a = particle->getDensity(),
                 P_a = particle->getPressure(),
                 c_a = particle->getSoundspeed();
    const point_t pos_a = particle->coordinates(),
                  vel_a = particle->getVelocity();


    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double rho_[n_nb],P_[n_nb],h_[n_nb],m_[n_nb],c_[n_nb],Pi_a_[n_nb];
    double vab_dot_DiWa_[n_nb];
    point_t pos_[n_nb], vel_[n_nb];

    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      rho_[b] = nb->getDensity();
      P_[b]   = nb->getPressure();
      pos_[b] = nb->coordinates();
      vel_[b] = nb->getVelocity();
      c_[b]   = nb->getSoundspeed();
      h_[b]   = nb->radius();
      m_[b]   = nb->mass() * (pos_[b]!=pos_a);
    }

    // precompute viscosity and kernel gradients
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double mu_ab = mu(h_a,h_[b],vel_a,vel_[b],pos_a,pos_[b]);
      Pi_a_[b] = artificial_viscosity(rho_a,rho_[b],c_a,c_[b],mu_ab);
      point_t        pos_ab = pos_a - pos_[b];
      space_vector_t vel_ab = point_to_vector(vel_a - vel_[b]);
      space_vector_t DiWab  = point_to_vector (
          kernel_gradient<kernel_type,gdimension>(pos_ab,.5*(h_a+h_[b])));
      vab_dot_DiWa_[b] = dot(vel_ab, DiWab);
    }

    // final answer
    double dudt_pressure, dudt_visc;
    dudt_visc = dudt_pressure = 0.0;
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      dudt_pressure += m_[b]*vab_dot_DiWa_[b];
      dudt_visc     += m_[b]*vab_dot_DiWa_[b]*Pi_a_[b];
    }
    double dudt = P_a/(rho_a*rho_a)*dudt_pressure + .5*dudt_visc;
    particle->setDudt(dudt);
  } // compute_dudt


  /**
   * @brief      Calculates the dedt, time derivative of either
   *             thermokinetic (internal + kinetic) or total
   *             (internal + kinetic + potential) energy.
   * See e.g. Rosswog (2009) "Astrophysical SPH" eq. (34)
   *
   *   de_a              (P_a*v_b   P_b*v_a   v_a + v_b       )
   *   ---- = -sum_b m_b ( -----  +  -----  + --------- Pi_ab ) . D_i Wab
   *    dt               (rho_a^2   rho_b^2       2           )
   *
   * @param      srch  The source's body holder
   * @param      nbsh  The neighbors' body holders
   */
  void compute_dedt(
      body* particle,
      std::vector<body*>& nbs)
  {


    using namespace viscosity;
    using namespace kernels;

    // this particle (index 'a')
    const double h_a = particle->radius(),
               rho_a = particle->getDensity(),
                 P_a = particle->getPressure(),
                 c_a = particle->getSoundspeed();
    const point_t pos_a = particle->coordinates(),
                  vel_a = particle->getVelocity();

    // neighbor particles (index 'b')
    const int n_nb = nbs.size();
    double rho_[n_nb],P_[n_nb],h_[n_nb],m_[n_nb],c_[n_nb],Pi_a_[n_nb];
    double va_dot_DiWa_[n_nb], vb_dot_DiWa_[n_nb];
    point_t pos_[n_nb], vel_[n_nb];

    for(int b = 0; b < n_nb; ++b) {
      const body * const nb = nbs[b];
      rho_[b] = nb->getDensity();
      P_[b]   = nb->getPressure();
      pos_[b] = nb->coordinates();
      vel_[b] = nb->getVelocity();
      c_[b]   = nb->getSoundspeed();
      h_[b]   = nb->radius();
      m_[b]   = nb->mass() * (pos_[b]!=pos_a);
    }

    // precompute viscosity and kernel gradients
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      double mu_ab = mu(h_a,h_[b],vel_a,vel_[b],pos_a,pos_[b]);
      Pi_a_[b] = artificial_viscosity(rho_a,rho_[b],c_a,c_[b],mu_ab);
      point_t        pos_ab = pos_a - pos_[b];
      space_vector_t DiWab  = point_to_vector (
          kernel_gradient<kernel_type,gdimension>(pos_ab,.5*(h_a+h_[b])));
      va_dot_DiWa_[b] = dot(point_to_vector(vel_a), DiWab);
      vb_dot_DiWa_[b] = dot(point_to_vector(vel_[b]), DiWab);
    }

    double dedt = 0;
    const double Prho2_a = P_a/(rho_a*rho_a);
    for(int b = 0 ; b < n_nb; ++b){ // Vectorized
      const double Prho2_b = P_[b]/(rho_[b]*rho_[b]);
      dedt -= m_[b]*( Prho2_a*vb_dot_DiWa_[b] + va_dot_DiWa_[b]*Prho2_b
               + .5*Pi_a_[b]*(vb_dot_DiWa_[b] + va_dot_DiWa_[b]));
    }
    particle->setDedt(dedt);
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
