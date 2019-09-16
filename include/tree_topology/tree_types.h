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
#include <iostream>
#include <math.h>
#include <vector>

#include "flecsi/geometry/point.h"

namespace flecsi {
namespace topology {

/*----------------------------------------------------------------------------*
 * class tree_branch
 * @brief Basic tree_branch implementation
 *----------------------------------------------------------------------------*/

//----------------------------------------------------------------------------//
//! \class tree_branch tree_types.h
//!
//! \brief tree_branch
//!
//! \tparam D Dimension
//! \tparam E Type for point (double)
//! \tparam KEY class of key used (hilbert, morton)
//----------------------------------------------------------------------------//
template <size_t D, typename E, class KEY> class tree_branch {
  using element_t = E;
  static constexpr size_t dimension = D;
  using key_t = KEY;
  using point_t = point_u<element_t, D>;

  //! Maximum number of children regarding the dimension
  static constexpr size_t num_children = 1 << dimension;

public:
  //! Describe the locality of a tree_branch
  enum b_locality : size_t { LOCAL = 0, EMPTY = 1, NONLOCAL = 2, SHARED = 3 };
  tree_branch() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    owner_ = rank;
  }
  tree_branch(const key_t &key) : key_(key) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    owner_ = rank;
  }
  tree_branch(const key_t &key, const point_t &coordinates,
              const element_t &mass, const point_t &bmin, const point_t &bmax,
              const b_locality &locality, const int &owner)
      : key_(key), coordinates_(coordinates), mass_(mass), bmin_(bmin),
        bmax_(bmax), locality_(locality), owner_(owner) {}

  ~tree_branch() { ents_.clear(); }

  // Getters
  key_t key() const { return key_; }
  point_t coordinates() { return coordinates_; };
  element_t mass() { return mass_; };
  point_t bmin() { return bmin_; };
  point_t bmax() { return bmax_; };
  uint64_t sub_entities() const { return sub_entities_; }
  b_locality locality() { return locality_; }
  int owner() { return owner_; };
  bool ghosts_local() { return ghosts_local_; };
  bool requested() { return requested_; };
  char bit_child() { return bit_child_; };

  // Setters
  void set_coordinates(const point_t &coordinates) {
    coordinates_ = coordinates;
  };
  void set_mass(const element_t &mass) { mass_ = mass; };
  void set_bmax(const point_t &bmax) { bmax_ = bmax; };
  void set_bmin(const point_t &bmin) { bmin_ = bmin; };
  void set_begin_tree_entities(const size_t &begin_tree_entities) {
    begin_tree_entities_ = begin_tree_entities;
  }
  void set_end_tree_entities(const size_t &end_tree_entities) {
    end_tree_entities_ = end_tree_entities;
  }
  size_t begin_tree_entities() { return begin_tree_entities_; }
  size_t end_tree_entities() { return end_tree_entities_; }
  void set_locality(b_locality locality) { locality_ = locality; }
  void set_sub_entities(uint64_t sub_entities) { sub_entities_ = sub_entities; }
  void set_leaf(bool leaf) { leaf_ = leaf; }
  void set_owner(int owner) { owner_ = owner; };
  void set_ghosts_local(const bool ghosts_local) {
    ghosts_local_ = ghosts_local;
  };
  void set_requested(bool requested) { requested_ = requested; };
  void set_bit_child(char bit_child) { bit_child_ = bit_child; };

  // Checker
  bool is_leaf() const { return leaf_; }
  bool is_valid() const { return true; }
  bool is_local() const {
    return locality_ == LOCAL || locality_ == EMPTY || locality_ == SHARED;
  }
  bool is_shared() const { return locality_ == SHARED; }

  //! Insert an entity in the branch vector
  void insert(const size_t &id) {
    assert(find(ents_.begin(), ents_.end(), id) == ents_.end());
    ents_.push_back(id);
  } // insert

  //! Number of entities in this branch
  int size() { return ents_.size(); }

  //! Remove a specific entity from the branch
  void remove(const size_t &id) {
    auto itr = find(ents_.begin(), ents_.end(), id);
    assert(itr != ents_.end());
    ents_.erase(itr);
  }
  // Remove a specific entity from the child bitmap
  void remove_bit(const int &bit) {
    assert(bit_child_ & (1 << bit));
    bit_child_ ^= 1 << bit;
    assert(!(bit_child_ & (1 << bit)));
  }

  auto begin() { return ents_.begin(); }
  auto end() { return ents_.end(); }
  auto clear() {
    ents_.clear();
    requested_ = false;
    ghosts_local_ = false;
  }

  //! Add a specific child in the child bitmap
  void add_bit_child(int i) {
    assert(!(bit_child_ & 1 << i));
    bit_child_ |= 1 << i;
  };
  //! Check if this branch have a specific child in the bitset
  bool as_child(int i) { return bit_child_ & 1 << i; };

protected:
  void set_key_(key_t key) { key_ = key; }

  friend std::ostream &operator<<(std::ostream &os, const tree_branch &b) {
    // TODO change regarding to dimension
    os << std::setprecision(10) << "Branch: coord: " << b.coordinates_;
    os << " mass: " << b.mass_ << " loc: " << b.locality_;
    os << " key: " << b.key_ << " sub_entities: " << b.sub_entities_;
    os << " owner: " << b.owner_
       << " bit_child: " << std::bitset<8>(b.bit_child_);
    return os;
  }

  key_t key_;                 // Key of this branch
  uint64_t sub_entities_ = 0; // Subentities in this subtree
  bool leaf_ = true;
  b_locality locality_ = EMPTY;
  point_t bmin_;
  point_t bmax_;
  int owner_;
  point_t coordinates_;
  element_t mass_;
  std::vector<size_t> ents_;
  bool ghosts_local_ = true;
  bool requested_ = false;
  char bit_child_ = 0;
  size_t begin_tree_entities_;
  size_t end_tree_entities_;
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
template <size_t D, typename E, class KEY> class entity {
  using element_t = E;
  static constexpr size_t dimension = D;
  using point_t = point_u<element_t, dimension>;
  using key_t = KEY;

public:
  entity(){};
  entity(const point_t &coordinates, const element_t &mass, const size_t &id,
         const element_t &radius, const key_t &key)
      : coordinates_(coordinates), mass_(mass), id_(id), radius_(radius),
        key_(key){};

  inline bool operator==(const entity &a) { return a.id_ == this->id_; }

  ~entity(){};
  // Getters
  point_t coordinates() const { return coordinates_; };
  element_t mass() const { return mass_; };
  key_t key() const { return key_; };
  element_t radius() const { return radius_; };
  size_t id() const { return id_; };
  int owner() { return owner_; }
  // Setters
  void set_coordinates(const point_t &coordinates) {
    coordinates_ = coordinates;
  };
  void set_mass(const element_t &mass) { mass_ = mass; };
  void set_radius(const element_t &radius) { radius_ = radius; };
  void set_key(const key_t &key) { key_ = key; };
  void set_id(const size_t &id) { id_ = id; };
  void set_owner(const int &owner) { owner_ = owner; };

  constexpr bool operator<(const entity &ent) const { return key_ <= ent.key_; }

  friend std::ostream &operator<<(std::ostream &os, const entity &b) {
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

/*----------------------------------------------------------------------------*
 * class tree_entity
 * @brief Basic tree_entity implementation
 *----------------------------------------------------------------------------*/

//----------------------------------------------------------------------------//
//! \class tree_entity tree_types.h
//!
//! \brief tree_entity
//!
//! \tparam D Dimension
//! \tparam E Type for point (double)
//! \tparam KEY class of key used (hilbert, morton)
//! \tparam ENT class of entities defined by the user
//----------------------------------------------------------------------------//
template <size_t D, typename E, class KEY, class ENT> class tree_entity {
public:
  using element_t = E;
  static constexpr size_t dimension = D;
  using point_t = point_u<element_t, dimension>;
  using key_t = KEY;

  using range_t = std::array<point_t, 2>;
  using entity_t = ENT;

protected:
  enum e_locality_ { LOCAL = 0, NONLOCAL = 1, SHARED = 2, EXCL = 3, GHOST = 4 };

public:
  tree_entity(const key_t &key, const point_t &coordinates,
              entity_t *entity_ptr, const int64_t owner, const element_t &mass,
              const size_t &id, const element_t &radius)
      : key_(key), coordinates_(coordinates), entity_ptr_(entity_ptr),
        owner_(owner), mass_(mass), id_(id), radius_(radius) {
    locality_ = entity_ptr_ == nullptr ? NONLOCAL : EXCL;
    global_id_ = id;
  };

  tree_entity() : key_(key_t::null()), locality_(NONLOCAL) {
    owner_ = -1;
    global_id_ = {};
  }

  tree_entity(const key_t &key, const point_t &coordinates)
      : key_(key), coordinates_(coordinates), locality_(NONLOCAL){};

  // Getters
  key_t key() const { return key_; }
  size_t id() const { return id_; }
  size_t global_id() const { return global_id_; }
  e_locality_ locality() { return locality_; }
  int64_t owner() { return owner_; }
  point_t coordinates() const { return coordinates_; }
  element_t mass() { return mass_; };
  element_t radius() { return radius_; };
  entity_t *entity_ptr() { return entity_ptr_; }

  // Setters
  void set_global_id(size_t id) { global_id_ = id; }
  void set_locality(e_locality_ loc) { locality_ = loc; }
  void set_owner(int64_t owner) { owner_ = owner; }
  void set_coordinates(point_t &coordinates) { coordinates_ = coordinates; }
  void set_radius(element_t radius) { radius_ = radius; };
  void set_shared() { locality_ = SHARED; };
  void set_entity_ptr(entity_t *entity_ptr) { entity_ptr_ = entity_ptr; }
  void set_id_(size_t id) { id_ = id; }
  void set_global_id_(size_t id) { global_id_ = id; }
  void set_entity_key_(key_t bid) { key_ = bid; }

  // Checker
  bool is_valid() const { return key_ != key_t::null(); }
  bool is_local() const {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return owner_ == rank;
  }

  friend std::ostream &operator<<(std::ostream &os, const tree_entity &b) {
    os << std::setprecision(10);
    os << "Tree entity. Pos: " << b.coordinates_ << " Mass: " << b.mass_ << " ";
    if (b.locality_ == LOCAL || b.locality_ == EXCL || b.locality_ == SHARED)
      os << "LOCAL";
    else
      os << "NONLOCAL";
    os << " owner: " << b.owner_ << " id: " << b.id_;
    os << " key: "<<b.key_;
    return os;
  }

private:
  element_t mass_;
  element_t radius_;
  point_t coordinates_;
  key_t key_;
  size_t id_;
  size_t global_id_;
  e_locality_ locality_;
  int64_t owner_;
  entity_t *entity_ptr_;
};

} // namespace topology
} // namespace flecsi
