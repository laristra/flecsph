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
 * @file physics.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic physics implementation
 */

#pragma once

#include <vector>

// Basic elements for the physics
namespace physics {
double dt = 0.0;
double dt_saved = 0.0;
double totaltime = 0.0;
double totaltime_next = 0.0;
double t_h5data_output = 0.0;
double t_screen_output = 0.0;
double t_scalar_output = 0.0;
int64_t iteration = 0;
} // namespace physics

#include "eforce.h"
#include "kernels.h"
#include "params.h"
#include "tree.h"
#include "user.h"
#include "utils.h"

#include "boundary.h"
#include "eos.h"
#include "integration.h"
#include "viscosity.h"
#include "tensor.h"
#include "fmm.h"

namespace physics {
using namespace param;

void
compute_cofm(node * cofm, std::vector<body *> ents, std::vector<node *> nodes) {
  // Then compute the CoFM
  point_t coordinates = point_t{};
  double radius = 0; // bmax
  double mass = 0;
  size_t sub_entities = 0;
  double lap = 0;
  point_t bmin, bmax;
  for(int i = 0; i < gdimension; ++i) {
    bmax[i] = -DBL_MAX;
    bmin[i] = DBL_MAX;
  }
  // Compute the center of mass and mass
  for(int i = 0; i < ents.size(); ++i) {
    body * ent = ents[i];
    // This correspond to a body
    coordinates += ent->mass() * ent->coordinates();
    mass += ent->mass();
    ++sub_entities;
    for(int d = 0; d < gdimension; ++d) {
      bmin[d] = std::min(bmin[d], ent->coordinates()[d] - ent->radius() / 2.);
      bmax[d] = std::max(bmax[d], ent->coordinates()[d] + ent->radius() / 2.);
    } // for
  }
  for(int i = 0; i < nodes.size(); ++i) {
    // This correspond to another node
    node * c = nodes[i];
    coordinates += c->mass() * c->coordinates();
    mass += c->mass();
    sub_entities += c->sub_entities();
    for(int d = 0; d < gdimension; ++d) {
      bmin[d] = std::min(bmin[d], c->bmin()[d]);
      bmax[d] = std::max(bmax[d], c->bmax()[d]);
    } // for
  } // for
  assert(mass != 0.);
  // Compute the radius
  coordinates /= mass;
  for(int i = 0; i < ents.size(); ++i) {
    body * ent = ents[i];
    double dist = distance(coordinates, ent->coordinates());
    radius = std::max(radius, dist);
    lap = std::max(lap, dist + ent->radius());
  }
  for(int i = 0; i < nodes.size(); ++i) {
    node * c = nodes[i];
    double dist = distance(coordinates, c->coordinates());
    radius = std::max(radius, dist + c->radius());
    lap = std::max(lap, dist + c->radius() + c->lap());
  } // for
  // Register and quit this node
  cofm->set_coordinates(coordinates);
  cofm->set_radius(radius);
  cofm->set_mass(mass);
  cofm->set_sub_entities(sub_entities);
  cofm->set_lap(lap);
  cofm->set_bmin(bmin);
  cofm->set_bmax(bmax);

}

/**
 * @brief      Subtracts mechanical energy from total energy
 *             to recover internal energy
 * @param      srch  The source's body holder
 */
void
recover_internal_energy(body & particle) {
  const point_t pos = particle.coordinates(),
                vel = particle.getVelocity();
  const double etot = particle.getTotalenergy(),
               ekin = .5*flecsi::dot(vel, vel),
               epot = external_force::potential(pos);
  const double eint = etot - ekin - epot;
  if (not (eint > 0)) {
    std::cerr << "ERROR: internal energy non-positive:" << std::endl
              << "particle id: " << particle.id()      << std::endl
              << "total energy: " << etot              << std::endl
              << "kinetic energy: " << ekin            << std::endl
              << "internal energy: " << eint           << std::endl
              << "potential energy: " << epot          << std::endl
              << "particle position: " << pos          << std::endl;
    mpi_assert(false);
  }
  particle.setInternalenergy(eint);
} // recover_internal_energy

/**
 * @brief      Using current internal energy and dudt,
 *             recompute pressure and soundspeed half-timestep ahead
 *
 * @param      particle  The particle body
 */
void
recompute_pressure_soundspeed(body& particle) {
  const double uint = particle.getInternalenergy();
  const double dudt = particle.getDudt();
  particle.setInternalenergy(uint + 0.5*dt*dudt);
  eos::compute_pressure(particle);
  eos::compute_soundspeed(particle);
  particle.setInternalenergy(uint);
}

/**
 * @brief      Using current total energy and dedt,
 *             recompute pressure and soundspeed half-timestep ahead
 *
 * @param      particle  The particle body
 */
void
recompute_pressure_soundspeed_thermokinetic(body& particle) {
  const double etot = particle.getTotalenergy();
  const double dedt = particle.getDedt();
  recover_internal_energy(particle);
  const double uint = particle.getInternalenergy();
  const point_t & v_a = particle.getVelocity();
  const point_t & a_a = particle.getAcceleration();
  const double v_dot_a = flecsi::dot(v_a, a_a);
  particle.setInternalenergy(uint + 0.5*dt*(dedt - v_dot_a));
  eos::compute_pressure(particle);
  eos::compute_soundspeed(particle);
  particle.setInternalenergy(uint);
}

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
compute_density(body & particle, std::vector<body *> & nbs) {
  using namespace kernels;
  const double h_a = particle.radius();
  const point_t pos_a = particle.coordinates();
  const int n_nb = nbs.size();
  mpi_assert(n_nb > 0);

  double r_a_[n_nb], m_[n_nb], h_[n_nb];
  for(int b = 0; b < n_nb; ++b) {
    const body * const nb = nbs[b];
    m_[b] = nb->mass();
    h_[b] = nb->radius();
    point_t pos_b = nb->coordinates();
    r_a_[b] = flecsi::magnitude(pos_a - pos_b);
  }

  double rho_a = 0.0;
  for(int b = 0; b < n_nb; ++b) { // Vectorized
    double Wab = sph_kernel_function(r_a_[b], .5 * (h_a + h_[b]));
    rho_a += m_[b] * Wab;
  } // for
  if(not(rho_a > 0)) {
    std::cout << "Density of a particle is not a positive number: "
              << "rho = " << rho_a << std::endl;
    std::cout << "Failed particle id: " << particle.id() << std::endl;
    std::cerr << "particle position: " << particle.coordinates() << std::endl;
    std::cerr << "particle velocity: " << particle.getVelocity() << std::endl;
    std::cerr << "particle acceleration: "
              << particle.getAcceleration() + particle.getGAcceleration()
              << std::endl;
    std::cerr << "smoothing length:  " << particle.radius() << std::endl;
    assert(false);
  }
  particle.setDensity(rho_a);
} // compute_density

/**
 * @brief      Computes maximum signal speed for the given particle
 *
 * @param      particle  The particle body
 * @param      nbs       Vector of neighbor particles
 */
void
compute_signalspeed(body & particle, std::vector<body *> & nbs) {
  using namespace param;
  using namespace kernels;
  using namespace flecsi;
  // this particle (index 'a')
  const double c_a = particle.getSoundspeed();
  const point_t pos_a = particle.coordinates(),
                  v_a = particle.getVelocity();

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  double c_a_[n_nb];
  point_t pos_[n_nb], n_a_[n_nb], v_a_[n_nb];

  for(int b = 0; b < n_nb; ++b) {
    const body * const nb = nbs[b];
    const point_t pos_b  = nb->coordinates();
    n_a_[b] = (pos_a - pos_b)/distance(pos_a, pos_b);
    v_a_[b]   = v_a - nb->getVelocity();
    c_a_[b]   = std::max(c_a, nb->getSoundspeed());
  }

  double vsig = 0.0;
  for(int b = 0 ; b < n_nb; ++b){
    vsig = std::max(vsig, c_a_[b] - std::min(dot(v_a_[b],n_a_[b]),0.0));
  }

  particle.setSignalspeed(vsig);
} // compute_signalspeed

/**
 * @brief      Compute divergence of the velocity field at this particle
 *
 * @param      particle  The particle body
 * @param      nbs       Vector of neighbor particles
 */
void
compute_divv(body & particle, std::vector<body *> & nbs) {
  using namespace kernels;

  // compute the divergence
  double div_v = 0.0;
  const double h_a = particle.radius();
  const point_t & pos_a = particle.coordinates();
  const point_t &   v_a = particle.getVelocity();
  for(int b = 0 ; b < nbs.size(); ++b){
    const body * const nb = nbs[b];
    const point_t & pos_b = nb->coordinates();
    const double h_ab = .5*(h_a + nb->radius());
    // const double h_b = nb->radius(); // DEBUG
    const double  m_b = nb->mass();
    const point_t DiWab = sph_kernel_gradient(pos_a - pos_b,h_ab);
    //point_t DiWab = .5*(sph_kernel_gradient(pos_a - pos_b,h_a)   // DEBUG
    //                 +  sph_kernel_gradient(pos_a - pos_b,h_b));
    div_v += m_b*dot(v_a, DiWab);
  }
  div_v /= particle.getDensity();

  // compute the divergence derivative
  const double div_v_p = particle.getDivergenceV();
  particle.setDdivvdt((div_v - div_v_p)/physics::dt);
  particle.setDivergenceV(div_v);

}

/**
 * @brief      Compute the density, EOS and soundspeed in one place
 * to save on gathering the neighbors
 *
 * @param      particle  The particle body
 * @param      nbs       Vector of neighbor particles
 */
void
compute_density_pressure_soundspeed(body & particle,
  std::vector<body *> & nbs) {
  compute_density(particle, nbs);
  if (evolve_internal_energy and thermokinetic_formulation)
    recover_internal_energy(particle);
  eos::compute_pressure(particle);
  eos::compute_soundspeed(particle);
  compute_signalspeed(particle, nbs);
  if (sph_viscosity == visc_cullen)
    compute_divv(particle,nbs);
}

/**
 * @brief      Calculates total energy for every particle
 *             NOTE: total energy does not include grav. energy
 * @param      particle
 */
void
set_total_energy(body & particle) {
  const point_t pos = particle.coordinates(),
                vel = particle.getVelocity();
  const double eint = particle.getInternalenergy(),
               ekin = .5*flecsi::dot(vel, vel),
               epot = external_force::potential(pos);
  particle.setTotalenergy(ekin + eint + epot);
} // set_total_energy

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
compute_acceleration(body & particle, std::vector<body *> & nbs) {
  using namespace param;
  using namespace viscosity;
  using namespace kernels;

  // Reset the acceleration

  // this particle (index 'a')
  const double h_a = particle.radius(),
             rho_a = particle.getDensity(),
               P_a = particle.getPressure(),
               c_a = particle.getSoundspeed(),
           alpha_a = particle.getAlpha();
  const point_t pos_a = particle.coordinates(),
                v12_a = particle.getVelocityhalf();

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  double rho_[n_nb],P_[n_nb],h_[n_nb],m_[n_nb],c_[n_nb],Pi_a_[n_nb],alpha_[n_nb];
  point_t pos_[n_nb], v12_[n_nb], DiWa_[n_nb];

  for(int b = 0; b < n_nb; ++b) {
    const body * const nb = nbs[b];
    rho_[b] = nb->getDensity();
    P_[b]   = nb->getPressure();
    pos_[b] = nb->coordinates();
    v12_[b] = nb->getVelocityhalf();
    c_[b]   = nb->getSoundspeed();
    h_[b]   = nb->radius();
    m_[b]   = nb->mass() * (pos_[b]!=pos_a); // if same particle, m_b->0
    alpha_[b] = nb->getAlpha();
  }

  // precompute viscosity and kernel gradients
  for(int b = 0; b < n_nb; ++b) { // Vectorized
    const point_t v12_ab = v12_a - v12_[b];
    const point_t pos_ab = pos_a - pos_[b];
    const double h_ab = .5*(h_a + h_[b]);
    const double mu_ab = mu(h_ab, v12_ab, pos_ab),
              alpha_ab = .5*(alpha_a + alpha_[b]),
                rho_ab = .5*(rho_a + rho_[b]),
                  c_ab = .5*(c_a + c_[b]);
    Pi_a_[b] = sph_artificial_viscosity(alpha_ab, rho_ab, c_ab, mu_ab);
    DiWa_[b] = sph_kernel_gradient(pos_ab,h_ab);
    // DiWa_[b] = .5*(sph_kernel_gradient(pos_ab,h_a)   // DEBUG
    //             + sph_kernel_gradient(pos_ab,h_[b]));
  }

  // compute the final answer
  const double Prho2_a = P_a / (rho_a * rho_a);
  point_t acc_a = 0.0;
  for(int b = 0; b < n_nb; ++b) { // Vectorized
    const double Prho2_b = P_[b] / (rho_[b] * rho_[b]);
    acc_a += -m_[b] * (Prho2_a + Prho2_b + Pi_a_[b]) * DiWa_[b];
  }
  acc_a += external_force::acceleration(particle);
  particle.setAcceleration(acc_a);
  particle.setGAcceleration(0);
  particle.setGPotential(0);
} // compute_acceleration

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
void
compute_dudt(body & particle, std::vector<body *> & nbs) {
  // Do not change internal energy in relaxation phase
  if(iteration < relaxation_steps) {
    particle.setDudt(0.0);
    return;
  }
  using namespace viscosity;
  using namespace kernels;

  // this particle (index 'a')
  const double h_a = particle.radius(),
             rho_a = particle.getDensity(),
               P_a = particle.getPressure(),
               c_a = particle.getSoundspeed(),
           alpha_a = particle.getAlpha();
  const point_t pos_a = particle.coordinates(),
                vel_a = particle.getVelocity(),
                v12_a = particle.getVelocityhalf();

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  double rho_[n_nb],P_[n_nb],h_[n_nb],m_[n_nb],c_[n_nb],Pi_a_[n_nb],alpha_[n_nb];
  double vab_dot_DiWa_[n_nb];
  point_t pos_[n_nb], vel_[n_nb], v12_[n_nb];

  for(int b = 0; b < n_nb; ++b) {
    const body * const nb = nbs[b];
    rho_[b] = nb->getDensity();
    P_[b]   = nb->getPressure();
    pos_[b] = nb->coordinates();
    vel_[b] = nb->getVelocity();
    v12_[b] = nb->getVelocityhalf();
    c_[b]   = nb->getSoundspeed();
    h_[b]   = nb->radius();
    m_[b]   = nb->mass() * (pos_[b]!=pos_a);
    alpha_[b] = nb->getAlpha();
  }

  // precompute viscosity and kernel gradients
  for(int b = 0 ; b < n_nb; ++b){ // Vectorized
    point_t pos_ab = pos_a - pos_[b];
    point_t v12_ab = v12_a - v12_[b];
    point_t vel_ab = vel_a - vel_[b];
    double h_ab = .5*(h_a + h_[b]);
    const double mu_ab = mu(h_ab, v12_ab, pos_ab),
              alpha_ab = .5*(alpha_a + alpha_[b]),
                rho_ab = .5*(rho_a + rho_[b]),
                  c_ab = .5*(c_a + c_[b]);
    Pi_a_[b] = sph_artificial_viscosity(alpha_ab, rho_ab, c_ab, mu_ab);
    point_t DiWab  = sph_kernel_gradient(pos_ab,h_ab);
    // point_t DiWab = .5*(sph_kernel_gradient(pos_ab,h_a)  // DEBUG
    //                  + sph_kernel_gradient(pos_ab,h_[b]));
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
  particle.setDudt(dudt);

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
void
compute_dedt(body & particle, std::vector<body *> & nbs) {
  using namespace viscosity;
  using namespace kernels;

  // this particle (index 'a')
  const double h_a = particle.radius(),
             rho_a = particle.getDensity(),
               P_a = particle.getPressure(),
               c_a = particle.getSoundspeed(),
           alpha_a = particle.getAlpha();
  const point_t pos_a = particle.coordinates(),
                vel_a = particle.getVelocity(),
                v12_a = particle.getVelocityhalf(),
                 ga_a = particle.getGAcceleration();
  const double gv = dot(ga_a,vel_a);

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  double rho_[n_nb],P_[n_nb],h_[n_nb],m_[n_nb],c_[n_nb],Pi_a_[n_nb],alpha_[n_nb];
  double va_dot_DiWa_[n_nb], vb_dot_DiWa_[n_nb];
  point_t pos_[n_nb], vel_[n_nb], v12_[n_nb];

  for(int b = 0; b < n_nb; ++b) {
    const body * const nb = nbs[b];
    rho_[b]   = nb->getDensity();
    P_[b]     = nb->getPressure();
    pos_[b]   = nb->coordinates();
    vel_[b]   = nb->getVelocity();
    v12_[b]   = nb->getVelocityhalf();
    c_[b]     = nb->getSoundspeed();
    h_[b]     = nb->radius();
    m_[b]     = nb->mass() * (pos_[b]!=pos_a);
    alpha_[b] = nb->getAlpha();
  }

  // precompute viscosity and kernel gradients
  for(int b = 0 ; b < n_nb; ++b){ // Vectorized
    point_t pos_ab = pos_a - pos_[b];
    point_t v12_ab = v12_a - v12_[b];
    point_t vel_ab = vel_a - vel_[b];
    double h_ab = .5*(h_a + h_[b]);
    const double mu_ab = mu(h_ab, v12_ab, pos_ab),
              alpha_ab = .5*(alpha_a + alpha_[b]),
                rho_ab = .5*(rho_a + rho_[b]),
                  c_ab = .5*(c_a + c_[b]);
    Pi_a_[b] = sph_artificial_viscosity(alpha_ab, rho_ab, c_ab, mu_ab);
    point_t DiWab = sph_kernel_gradient(pos_ab,h_ab);
    // point_t DiWab = .5*(sph_kernel_gradient(pos_ab,h_a) // DEBUG
    //                  + sph_kernel_gradient(pos_ab,h_[b]));
    va_dot_DiWa_[b] = dot(vel_a, DiWab);
    vb_dot_DiWa_[b] = dot(vel_[b], DiWab);
  }

  double dedt = 0;
  const double Prho2_a = P_a/(rho_a*rho_a);
  for(int b = 0 ; b < n_nb; ++b){ // Vectorized
    const double Prho2_b = P_[b]/(rho_[b]*rho_[b]);
    dedt -= m_[b]*( Prho2_a*vb_dot_DiWa_[b] + va_dot_DiWa_[b]*Prho2_b
             + .5*Pi_a_[b]*(vb_dot_DiWa_[b] + va_dot_DiWa_[b]));
  }
  dedt += gv;
  particle.setDedt(dedt);

} // compute_dedt

/**
 * @brief      Adds energy dissipation rate due to artificial
 *             particle relaxation drag force
 * @param      particle
 */
void
add_drag_dedt(body & particle) {
  using namespace param;
  const point_t vel = particle.getVelocity();
  const point_t acc = external_force::acceleration_drag(vel);
  const double v_dot_a = flecsi::dot(vel, acc);
  double dedt = particle.getDedt();
  particle.setDedt(dedt + v_dot_a);
} // add_drag_dedt

/**
 * @brief      Adds energy dissipation rate due to artificial
 *             particle relaxation drag force - internal energy
 * @param      srch  The source's body holder
 */
void add_drag_dudt(body& source) {
  using namespace param;
  const point_t vel = source.getVelocity();
  const point_t acc = external_force::acceleration_drag(vel);
  const double v_dot_a = flecsi::dot(vel, acc);
  double dudt = source.getDudt();
  source.setDudt(dudt + v_dot_a);
} // add_drag_dudt

/**
 * @brief      Compute the timestep from acceleration and mu
 *
 * @param      particle
 */
void
compute_dt(body & particle) {
  const double tiny = 1e-24;
  const double mc   = 0.6; // constant in denominator for viscosity

  // particles separation around this particle
  const double dx = particle.radius()
                  / (sph_eta*kernels::kernel_width);

  // timestep based on particle velocity
  const point_t vel = particle.getVelocity();
  const double vn  = magnitude(vel);
  const double dt_v = dx/(vn + tiny);

  // timestep based on acceleration
  const double acc = magnitude(particle.getAcceleration()
                             + particle.getGAcceleration());
  const double dt_a = sqrt(dx/(acc + tiny));

  // timestep based on sound speed and viscosity
  const double cs_a = particle.getSoundspeed();
  const double vsig_a = particle.getSignalspeed();
  const double dt_c = dx/(tiny + vsig_a*(1 + mc*sph_viscosity_alpha));

  // minimum timestep
  double dtmin = timestep_cfl_factor * std::min(std::min(dt_v,dt_a), dt_c);

  // timestep based on positivity of internal energy
  if (evolve_internal_energy and thermokinetic_formulation) {
    const point_t pos = particle.coordinates(),
                  gra = particle.getGAcceleration();
    const double eint = particle.getInternalenergy(),
                 epot = external_force::potential(pos);
    double delta_epot, delta_egrv;
    if (not (eint > 0)) {
      std::cerr
          << "ERROR: internal energy non-positive "
          << "for particle " << particle.id() << std::endl
          << "particle position: " << pos   << std::endl
          << "particle velocity: " << vel   << std::endl
          << "particle acceleration: "
          << particle.getAcceleration() + particle.getGAcceleration()
          << "internal energy: "   << eint  << std::endl
          << "potential energy: "  << epot  << std::endl
          << "gravitationial potential: "
          << particle.getGPotential() << std::endl
          << "total energy: "<< particle.getTotalenergy() << std::endl
          << "smoothing length:  " << particle.radius()   << std::endl;
      assert(false);
    }

    int i;
    for(i=0; i<20; ++i) { // '20' hardcoded parameter (# or iterations)
      delta_epot =  external_force::potential(pos + dtmin*vel) - epot;
      delta_egrv = -dtmin*flecsi::dot(vel, gra);
      if(delta_epot + delta_egrv < eint*0.5) break; // '0.5' hardcoded
      dtmin *= 0.25;                                // '0.25' hardcoded
    }

    if (i>=20) {
      std::cerr << "ERROR: eint-based dt estimator loop did not converge "
                << "for particle " << particle.id() << std::endl;
      std::cerr << "particle position: " << pos << std::endl
                << "particle velocity: " << vel << std::endl
                << "particle acceleration: "
                << particle.getAcceleration() + particle.getGAcceleration()
                << std::endl;
      std::cerr << "smoothing length:  " << particle.radius()
                                         << std::endl;
      std::cerr << "dx: " << dx << std::endl;
      std::cerr << "dt_v = " << dt_v << ", vn = " << vn << std::endl;
      std::cerr << "dt_a = " << dt_a << ", acc = " << acc << std::endl;
      std::cerr << "dt_c = " << dt_c << ", cs_a = " << cs_a << std::endl;
      std::cerr << "internal energy: " << eint << std::endl;
      std::cerr << "potential energy: " << epot << std::endl;
      std::cerr << "gravitationial potential: "
                << particle.getGPotential() << std::endl;
      std::cerr << "total energy: "<< particle.getTotalenergy() << std::endl;
      dtmin = timestep_cfl_factor * std::min(std::min(dt_v,dt_a), dt_c);
      for(i=0; i<20; ++i) { // '20' hardcoded parameter (# or iterations)
        delta_epot =  external_force::potential(pos + dtmin*vel) - epot;
        delta_egrv = -dtmin*flecsi::dot(vel, gra);
        std::cerr << "dtmin[" << i << "] = " << dtmin
                  << ", epot = " << epot + delta_epot
                  << ", delta_egrv = " << delta_egrv << std::endl;
        if(delta_epot + delta_egrv < eint*0.5) break; // '0.5' hardcoded
        dtmin *= 0.25;                                // '0.25' hardcoded
      }
      assert(false);
    }
  }

  particle.setDt(dtmin);
}

/**
 * @brief      Adds dissipative drag to acceleration
 *             (used in particles relaxation step)
 * @param      particle
 */
void
add_drag_acceleration(body & particle) {
  using namespace param;
  point_t acc = particle.getAcceleration();
  const point_t vel = particle.getVelocity();
  acc += external_force::acceleration_drag(vel);
  particle.setAcceleration(acc);
} // add_drag_acceleration

/**
 * @brief      Short-range repulsion force
 *
 *     (dv_a)                             (  r_b - r_a       )
 *     (----)   += -gamma_repulsion  sum_b( ---------- * m_b )
 *     ( dt )_i                           (  |r_ab|^3        )
 *
 *             Artificial force to prevent particles from clumping
 *
 * @param      particle  The particle body
 * @param      nbs       Vector of neighbor particles
 */
void
add_short_range_repulsion(body & particle, std::vector<body *> & nbs) {
  using namespace param;

  // this particle (index 'a')
  const double h_a = particle.radius();
  const size_t id_a = particle.id();
  const point_t pos_a = particle.coordinates();
  point_t acc_a = particle.getAcceleration();

  // neighbor particles (index 'b')
  const int n_nb = nbs.size();
  double h_b, m_b;
  point_t pos_b;
  point_t acc_r = 0.0;

  for(int b = 0; b < nbs.size(); ++b) {
    const body * const nb = nbs[b];
    if(nb->id() == id_a)
      continue;
    h_b = nb->radius();
    pos_b = nb->coordinates();
    double h_ab = .5 * (h_a + h_b);
    double r_ab = flecsi::magnitude(pos_a - pos_b);
    if(r_ab > h_ab * relaxation_repulsion_radius)
      continue;
    m_b = nb->mass();
    acc_r += m_b * (pos_a - pos_b) / (r_ab * r_ab * r_ab);
  }
  acc_r *= relaxation_repulsion_gamma;
  particle.setAcceleration(acc_a + acc_r);
} // add_short_range_repulsion

/**
 * @brief      Reduce adaptive timestep and set its value
 *
 * @param      bodies   Set of bodies
 */
void
set_adaptive_timestep(std::vector<body> & bodies) {
  double dtmin = 1e24; // some ludicrous number

  for(size_t i = 0; i < bodies.size(); ++i) {
    dtmin = std::min(dtmin, bodies[i].getDt());
  }

  mpi_utils::reduce_min(dtmin);

  if(dtmin < dt)
    dt = std::min(dtmin, dt / 2.0);

  if(dtmin > 2.0 * dt)
    dt = dt * 2.0;

  double dt_orig = dt;
  totaltime_next = totaltime + dt;

  if(out_screen_dt > 0) {
    // if output to screen by time:
    // match the next screen output time
    if (totaltime == t_screen_output) {
      if (dt_saved > 0)
        dt = dt_saved;
      t_screen_output += out_screen_dt;
    }
    if (totaltime + dt > t_screen_output - 0.1*out_screen_dt) {
      dt_saved = dt;
      totaltime_next = t_screen_output;
    }
  }

  if(out_scalar_dt > 0) {
    // if output scalar by time:
    // match the next scalar output time
    if (totaltime == t_scalar_output) {
      if (dt_saved > 0)
        dt = dt_saved;
      t_scalar_output += out_scalar_dt;
    }
    if (totaltime + dt > t_scalar_output - 0.1*out_scalar_dt) {
      dt_saved = dt;
      totaltime_next = std::min(totaltime_next, t_scalar_output);
    }
  }

  if(out_h5data_dt > 0) {
    // if output h5data by time:
    // match the next h5data output time
    if (totaltime == t_h5data_output) {
      if (dt_saved > 0)
        dt = dt_saved;
      t_h5data_output += out_h5data_dt;
    }
    if (totaltime + dt > t_h5data_output - 0.1*out_h5data_dt) {
      dt_saved = dt;
      totaltime_next = std::min(totaltime_next, t_h5data_output);
    }
  }
  dt = totaltime_next - totaltime;

}

void
compute_smoothinglength(std::vector<body> & bodies) {
  if constexpr (gdimension == 1) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      double m_b = bodies[i].mass();
      double rho_b = bodies[i].getDensity();
      bodies[i].set_radius(m_b / rho_b * sph_eta * kernels::kernel_width);
    }
  }
  if constexpr (gdimension == 2) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      double m_b = bodies[i].mass();
      double rho_b = bodies[i].getDensity();
      bodies[i].set_radius(sqrt(m_b / rho_b) * sph_eta * kernels::kernel_width);
    }
  }
  if constexpr (gdimension == 3) {
    for(size_t i = 0; i < bodies.size(); ++i) {
      double m_b = bodies[i].mass();
      double rho_b = bodies[i].getDensity();
      bodies[i].set_radius(cbrt(m_b / rho_b) * sph_eta * kernels::kernel_width);
    }
  } // if gdimension
}

/**
 * @brief  Advance time
 */
void
advance_time() {
  iteration++;
  if (adaptive_timestep)
    totaltime = totaltime_next;
  else
    totaltime = totaltime + dt;
}

/**
 * @brief  Termination criteria
 */
bool
termination_criteria() {
    if (final_iteration > 0 && iteration > final_iteration)
       return true;
    if (final_time > 0.0 && totaltime > final_time)
       return true;
    return false;
}

/**
 * @brief update smoothing length for particles (Rosswog'09, eq.51)
 *
 * ha = eta/N \sum_b pow(m_b / rho_b,1/dimension)
 */
void
compute_average_smoothinglength(std::vector<body> & bodies,
  int64_t nparticles) {
  compute_smoothinglength(bodies);
  // Compute the total
  double total = 0.;

  for(size_t i = 0; i < bodies.size(); ++i) {
    total += bodies[i].radius();
  }

  // Add up with all the processes
  MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Compute the new smoothing length
  double new_h = 1. / (double)nparticles * total;
  for(size_t i = 0; i < bodies.size(); ++i) {
    bodies[i].set_radius(new_h);
  }
}

/**
 * @brief      Checks all fields of the particle for NaNs
 * @param      particle  The particle to be checked
 */
#define NANCHECK_POINT_T(vfield) \
  { auto vfield = particle.vfield(); \
  for (int d = 0; d < gdimension; ++d) { \
    if (vfield[d] != vfield[d]) { \
      log_one(error) \
          << "particle[" << id << "]: NaN in " #vfield << std::endl; \
      passed = false; }}}
#define NANCHECK_DOUBLE(dfield) \
  { auto dfield = particle.dfield(); \
    if (dfield != dfield) { \
      log_one(error) \
          << "particle[" << id << "]: NaN in " #dfield << std::endl; \
      passed = false; }}

void
check_nans(body & particle) {
  auto id = particle.id();
  bool passed = true;
  if (id != id) {
    log_one(error) << "particle id is NaN: " << id << std::endl;
    passed = false;
  }
  NANCHECK_POINT_T(coordinates)
  NANCHECK_POINT_T(getVelocity)
  NANCHECK_POINT_T(getVelocityhalf)
  NANCHECK_POINT_T(getAcceleration)
  NANCHECK_POINT_T(getGAcceleration)
  NANCHECK_DOUBLE(mass)
  NANCHECK_DOUBLE(getGPotential)
  NANCHECK_DOUBLE(getDensity)
  NANCHECK_DOUBLE(getPressure)
  NANCHECK_DOUBLE(getEntropy)
  NANCHECK_DOUBLE(getTotalenergy)
  NANCHECK_DOUBLE(getDedt)
  NANCHECK_DOUBLE(getDudt)
  NANCHECK_DOUBLE(getAdiabatic)
  NANCHECK_DOUBLE(getSignalspeed)
  if (evolve_internal_energy) {
    NANCHECK_DOUBLE(getInternalenergy)
  }
  assert (passed);
} // check_nans
#undef NANCHECK_DOUBLE
#undef NANCHECK_POINT_T

/**
 * @brief      Stops simulation if any negativity is detected
 * @param      particle  The particle to be checked
 */
void
check_negativity(body & particle) {
  auto id  = particle.id();
  auto rho = particle.getDensity();
  auto P   = particle.getPressure();
  auto u   = particle.getInternalenergy();
  bool passed = true;
  if (rho < 0) {
    log_one(error)
        << "particle[" << id << "]: negative density = "
        << rho << std::endl;
    passed = false;
  }
  if (P < 0) {
    log_one(error)
        << "particle[" << id << "]: negative pressure = "
        << P << std::endl;
    passed = false;
  }
  if (param::evolve_internal_energy and u < 0) {
    log_one(error)
        << "particle[" << id << "]: negative internal energy = "
        << u << std::endl;
    passed = false;
  }
  assert (passed);
} // check_negativity

}; // namespace physics

