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
 * @file tree_types.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Representation of a body, a particle for our SPH implementation
 */

#pragma once

#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>

#include <mutex>

#include "space_vector.h"

namespace flecsi {
namespace topology {

enum type : char { NODE = 0, ENTITY = 1 };

/*----------------------------------------------------------------------------*
 * class cofm_u
 * @brief Basic center of mass implementation
 *----------------------------------------------------------------------------*/

template<size_t D, typename E, class KEY>
class cofm_u
{
  using element_t = E;
  using point_t = space_vector_u<element_t, D>;
  using key_t = KEY;

public:
  cofm_u() {
    coordinates_ = point_t{};
    mass_ = 0.;
    sub_entities_ = 0;
    radius_ = 0.;
  };

  cofm_u(const key_t & key) : key_(key) {
    coordinates_ = point_t{};
    mass_ = 0.;
    sub_entities_ = 0;
    radius_ = 0.;
    bmin_ = point_t{};
    bmax_ = point_t{};
  };

  cofm_u(const cofm_u & c) {
    coordinates_ = c.coordinates();
    mass_ = c.mass();
    radius_ = c.radius();
    sub_entities_ = c.sub_entities();
    lap_ = c.lap();
    key_ = c.key();
    bmin_ = c.bmin();
    bmax_ = c.bmax();
  }

  point_t coordinates() const {
    return coordinates_;
  }
  element_t mass() const {
    return mass_;
  }
  element_t radius() const {
    return radius_;
  }
  int sub_entities() const {
    return sub_entities_;
  }
  element_t lap() const {
    return lap_;
  }
  key_t key() const {
    return key_;
  }
  point_t bmin() const {
    return bmin_;
  }
  point_t bmax() const {
    return bmax_;
  }

  void set_coordinates(const point_t & coordinates) {
    coordinates_ = coordinates;
  }
  void set_mass(const element_t & mass) {
    mass_ = mass;
  }
  void set_radius(const element_t & radius) {
    radius_ = radius;
  }
  void set_sub_entities(const int & sub_entities) {
    sub_entities_ = sub_entities;
  }
  void set_lap(const element_t & lap) {
    lap_ = lap;
  }
  void set_bmin(const point_t & bmin) {
    bmin_ = bmin;
  }
  void set_bmax(const point_t & bmax) {
    bmax_ = bmax;
  }

protected:
  point_t coordinates_;
  element_t mass_;
  element_t radius_;
  point_t bmin_, bmax_;
  int sub_entities_;
  element_t lap_;
  key_t key_;
}; // class cofm

/**
 * @brief Class hcell, a cell in the hashtable
 * that represents the tree topology
 **/
template<size_t D, class KEY, class NODE, class ENTITY>
class hcell
{
  static constexpr int dimension = D;
  static constexpr int nchildren_ = 1 << dimension;
  using key_t = KEY;

  enum type_displ : int {
    CHILD_DISPL = 0,
    LOCALITY_DISPL = 1 << dimension,
    REQUESTED_DISPL = (1 << dimension) + 2,
    NCHILD_RECV_DISPL = (1 << dimension) + 3
  };
  enum type_mask : int {
    CHILD_MASK = 0b11111111,
    LOCALITY_MASK = 0b11 << LOCALITY_DISPL,
    REQUESTED_MASK = 0b1 << REQUESTED_DISPL,
    NCHILD_RECV_MASK = 0b1111 << NCHILD_RECV_DISPL
  };
  enum type_locality : int { LOCAL = 0, NONLOCAL = 1, SHARED = 2 };

public:
  hcell(const key_t & key) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    owner_ = rank_;
    key_ = key;
    node_idx_ = -1;
    entity_idx_ = -1;
    type_ = 0;
  }

  hcell(const key_t & key, const int entity_idx) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    owner_ = rank_;
    key_ = key;
    node_idx_ = -1;
    entity_idx_ = entity_idx;
    type_ = 0;
  }

  bool get_child(const int & c) const {
    return type_ & (1 << c);
  }
  void add_child(const int & c) {
    type_ = type_ | (1 << c);
  }
  int nchildren() const {
    int nchild = 0;
    for(int i = 0; i < nchildren_; ++i)
      nchild += get_child(i);
    return nchild;
  }
  void set_node_idx(const int node_idx) {
    node_idx_ = node_idx;
    assert(entity_idx_ == -1);
  }
  void set_entity_idx(const int entity_idx) {
    entity_idx_ = entity_idx;
    assert(node_idx_ == -1);
  }
  void set_shared() {
    type_ &= ~LOCALITY_MASK;
    type_ |= SHARED << LOCALITY_DISPL;
  }
  void set_requested() {
    type_ &= ~REQUESTED_MASK;
    type_ |= REQUESTED_MASK;
  }
  void unset_requested() {
    type_ &= ~REQUESTED_MASK;
  }
  /*
   * Number of children the cell is expected to receive
   */
  int nchildren_to_receive() const {
    return (type_ >> NCHILD_RECV_DISPL) & 0b1111;
  }
  void set_nchildren_to_receive(const int n) {
    type_ &= ~NCHILD_RECV_MASK;
    type_ |= (n << NCHILD_RECV_DISPL);
  }

  void set_owner(const int & owner) {
    owner_ = owner;
  }

  bool iam_owner() const {
    return owner_ == rank_;
  }

  bool is_shared() const {
    return ((type_ & LOCALITY_MASK) >> LOCALITY_DISPL) == SHARED;
  }

  bool requested() {
    return (type_ & REQUESTED_MASK);
  }

  bool is_empty_node() const {
    return is_node() && (!has_child() || (nchildren_to_receive() > nchildren()));
  }
  bool has_child() const {
    return type_ & (1 << (1 << dimension)) - 1;
  }

  int node_idx() const {
    return node_idx_;
  }
  int entity_idx() const {
    return entity_idx_;
  }
  unsigned int type() const {
    return type_;
  }
  key_t key() const {
    return key_;
  }
  int owner() const {
    return owner_;
  }

  bool is_node() const {
    assert(node_idx_ != -1 || entity_idx_ != -1);
    return node_idx_ != -1;
  }
  bool is_entity() const {
    return !is_node();
  }

  bool is_unset() const {
    return node_idx_ == -1 && entity_idx_ == -1;
  }

private:
  KEY key_;
  int node_idx_ = -1;
  int entity_idx_ = -1;
  int owner_;
  unsigned int type_ = 0;
  int rank_;
};

/*----------------------------------------------------------------------------*
 * class entity
 * @brief Basic entity implementation
 *----------------------------------------------------------------------------*/

//----------------------------------------------------------------------------//
//! \class entity tree_types.h
//!
//! \brief entity
//!
//! \tparam D Dimension
//! \tparam E Type for point (double)
//! \tparam KEY class of key used (hilbert, morton)
//----------------------------------------------------------------------------//
template<size_t D, typename E, class KEY>
class entity
{
  using element_t = E;
  static constexpr size_t dimension = D;
  using point_t = space_vector_u<element_t, dimension>;
  using key_t = KEY;

public:
  entity(){};
  entity(const point_t & coordinates,
    const element_t & mass,
    const size_t & id,
    const element_t & radius,
    const key_t & key)
    : coordinates_(coordinates), mass_(mass), id_(id), radius_(radius),
      key_(key){};

  inline bool operator==(const entity & a) {
    return a.id_ == this->id_;
  }

  ~entity(){};
  // Getters
  point_t coordinates() const {
    return coordinates_;
  };
  element_t mass() const {
    return mass_;
  };
  key_t key() const {
    return key_;
  };
  element_t radius() const {
    return radius_;
  };
  size_t id() const {
    return id_;
  };
  int owner() {
    return owner_;
  }
  // Setters
  void set_coordinates(const point_t & coordinates) {
    coordinates_ = coordinates;
  };
  void set_mass(const element_t & mass) {
    mass_ = mass;
  };
  void set_radius(const element_t & radius) {
    radius_ = radius;
  };
  void set_key(const key_t & key) {
    key_ = key;
  };
  void set_id(const size_t & id) {
    id_ = id;
  };
  void set_owner(const int & owner) {
    owner_ = owner;
  };

  constexpr bool operator<(const entity & ent) const {
    return key_ <= ent.key_;
  }

  friend std::ostream & operator<<(std::ostream & os, const entity & b) {
    // TODO change regarding to dimension
    os << std::setprecision(10) << "Particle: coord: " << b.coordinates_;
    os << " mass: " << b.mass_ << " h: " << b.radius_ << " id: " << b.id_;
    os << "key: " << b.key_ << " owner: " << b.owner_;
    return os;
  }

protected:
  point_t coordinates_;
  element_t mass_;
  size_t id_;
  element_t radius_;
  key_t key_;
  int owner_;
}; // class entity

} // namespace topology
} // namespace flecsi
