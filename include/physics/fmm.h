/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
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
 * @file fmm.h
 * @brief Functions used in the FMM computation
 */

#pragma once

#include "params.h"
#include "tree.h"

namespace fmm {
using namespace param;
double gc = gravitational_constant;

/*
 * @brief Compute gravitation interaction between two points
 *        Returns the resulting gravitational acceleration
 */
inline point_t 
gravitation_p2p(double & gpot,
  const point_t & local_coordinates,
  const point_t & dist_coordinates,
  const double & sm) {
  double dist = flecsi::distance(local_coordinates,dist_coordinates);
  gpot += -gc*sm/dist;
  point_t res =-gc*sm/(dist*dist*dist)*(local_coordinates - dist_coordinates);
  return res;
}

/*
 * @brief Compute the gravitation interaction between point and cell
 */
inline void 
gravitation_fc(double & pc,
  point_t & fc,
  const point_t & local_coordinates,
  const node * source) {
  const point_t & dist_coordinates = source->coordinates();
  const double M = source->mass();
  double d = flecsi::distance(local_coordinates,dist_coordinates);
  double d3 = d*d*d;
  point_t r = local_coordinates - dist_coordinates;

  pc += -gc*M/d;
  for(int m = 0; m < gdimension; ++m) {
    fc[m] += -gc*M*r[m]/d3; // Monopole
  }
}

/*
 * @brief Compute the gravitation interaction between point and node
 */
inline void
gravitation_fc(double & pc,
  point_t & fc,
  const point_t & local_coordinates,
  const body * source) {
  const point_t & dist_coordinates = source->coordinates();
  const double M = source->mass();
  double d = flecsi::distance(local_coordinates,dist_coordinates);
  double d3 = d*d*d;
  point_t r = local_coordinates-dist_coordinates;
  pc += -gc*M/d;
  for(int m = 0; m < gdimension; ++m){
    fc[m] += -gc*M*r[m]/d3; // Monopole
  }
}

/*
 * @brief Taylor expansion up to first order using gravity 
 *        at the cell center of mass
 */
void 
interaction_c2p(body * sink, const node * source) {
  const double & pc  = source->pc();
  const point_t & fc = source->fc();
  point_t cofm_coordinates = source->coordinates();

  point_t part_coordinates = sink->coordinates();
  point_t r = part_coordinates - cofm_coordinates;
  point_t grav = fc;
  double pot = pc;

  for(int i = 0 ; i < gdimension; ++i){
    pot += -r[i]*fc[i];
  }
  sink->setGPotential(sink->getGPotential()+pot);
  sink->setGAcceleration(grav+sink->getGAcceleration());
}


/**
 * @brief For all sub_entities, add gravitational force and potential of the
 *        node nd, using Taylor expansion coefficients stored in the node.
 */
void
fmm_c2p(const node* nd, std::vector<body *> & sub_entities) {
  for (int k = 0; k < sub_entities.size(); ++k) {
    interaction_c2p(sub_entities[k], nd);
  }
}

/**
 * @brief Particle-particle interactions between 'sources' and 'sinks'
 */
void 
fmm_p2p(std::vector<body *> & sinks,
  const node * node_sources,
  const std::vector<body *> & particle_sources) {
  for (int i=0; i<sinks.size(); ++i) {
    body *p = sinks[i];
    double pc = p->getGPotential();
    point_t acc = p->getGAcceleration();
    if(node_sources != nullptr)
    gravitation_fc(pc, acc, p->coordinates(), node_sources);
    for (int k = 0; k < particle_sources.size(); ++k) {
      body *q = particle_sources[k];
      if (q->id() == p->id())
        continue;
      acc+= gravitation_p2p(pc,p->coordinates(),q->coordinates(),q->mass());
    } // for
    p->setGPotential(pc);
    p->setGAcceleration(acc);
  }
}

/**
 * @brief node-node interaction: update Taylor expansion coefficients
 */
void
taylor_c2c(node * sink, const node * source) {
  gravitation_fc(sink->pc(), sink->fc(), sink->coordinates(), source);
}

/**
 * @brief node<-particle interaction: update Taylor expansion coefficients
 */
void 
taylor_p2c(node * sink, const body * source) {
  gravitation_fc(sink->pc(), sink->fc(), sink->coordinates(), source);
}

} // namespace fmm
