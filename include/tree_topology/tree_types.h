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

#include <mutex>

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
  static constexpr uint dimension = D;
  using key_t = KEY;
  using point_t = point_u<element_t, D>;

  //! Maximum number of children regarding the dimension
  static constexpr uint num_children = 1 << dimension;

public:
  //! Describe the locality of a tree_branch
  enum b_locality : uint { EMPTY = 0, LOCAL = 1, NONLOCAL = 2, SHARED = 3 };

  enum bits : uint {CHILD = 0, LEAF = 8, LOC = 9, GHS_LOC = 11, RQST = 12};

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
        bmax_(bmax), owner_(owner) {
    set_locality(locality); 
  }

  ~tree_branch() { ents_.clear(); }

  // Getters
  key_t key() const { return key_; }
  point_t coordinates() const { return coordinates_; };
  element_t mass() const { return mass_; };
  point_t bmin() const { return bmin_; };
  point_t bmax() const { return bmax_; };
  uint64_t sub_entities() const { return sub_entities_; }
  int owner() const { return owner_; };

  b_locality locality() const {  
    return static_cast<b_locality>((type_ & (uint(3) << bits::LOC)) >> bits::LOC);
  }
  bool ghosts_local() const { return type_ & (uint(1) << bits::GHS_LOC);}
  bool requested() const { return type_ & (uint(1) << bits::RQST);}
  char bit_child() const { return static_cast<char>(type_ & uint(255));}
  bool is_leaf() const { 
    bool res = type_ & (uint(1) << bits::LEAF);
    return type_ & (uint(1) << bits::LEAF); 
  }
  bool is_local() const {
    b_locality loc = locality(); 
    return loc == LOCAL || loc == EMPTY || loc == SHARED;
  }
  bool is_shared() const {
    b_locality loc = locality(); 
    return loc == SHARED; 
  }
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
  size_t begin_tree_entities()const { return begin_tree_entities_; }
  size_t end_tree_entities() const { return end_tree_entities_; }
  void set_sub_entities(const uint64_t& sub_entities) { sub_entities_ = sub_entities; }
  void set_owner(const int& owner) { owner_ = owner; };

  void set_locality(const b_locality& locality) { 
    uint value = locality; 
    assert(value >= 0 || value <= 3); 
    constexpr uint mask_loc = ~(uint(3) << bits::LOC);  
    type_ &= mask_loc;
    type_ |= value << bits::LOC; 
  }
  void set_leaf(const bool leaf) { 
    uint value = leaf; 
    assert(value == 0 || value == 1); 
    constexpr uint mask_leaf = ~(uint(1) << bits::LEAF); 
    type_ &= mask_leaf; 
    type_ |= value << bits::LEAF; 
  }
  void set_ghosts_local(const bool ghosts_local) {
    uint value = ghosts_local; 
    assert(value == 0 || value == 1); 
    constexpr uint mask_gl = ~(uint(1) << bits::GHS_LOC); 
    type_ &= mask_gl; 
    type_ |= value << bits::GHS_LOC; 
  }
  void set_requested(const bool requested) {
    uint value = requested; 
    assert(value == 0 || value == 1); 
    constexpr uint mask_rq = ~(uint(1) << bits::RQST); 
    type_ &= mask_rq; 
    type_ |= value << bits::RQST; 
  }
  void set_bit_child(const uint bit_child) { 
    constexpr uint mask_ch = ~uint(255); 
    type_ &= mask_ch;
    type_ |= bit_child; 
  }
  //! Add a specific child in the child bitmap
  void add_bit_child(const uint i) {
    assert(!(type_ & uint(1) << i));
    type_ |= uint(1) << i;
  }
  //! Check if this branch have a specific child in the bitset
  bool as_child(const uint i) { 
    bool res = (type_ & (uint(1)<<i)); 
    return type_ & uint(1) << i; 
  };
  // Remove a specific entity from the child bitmap
  void remove_bit(const uint& bit) {
    assert(type_ & (1 << bit));
    type_ ^= uint(1) << bit;
    assert(!(type_ & (uint(1) << bit)));
  }

  //! Insert an entity in the branch vector
  void insert(const size_t &id) {
    assert(find(ents_.begin(), ents_.end(), id) == ents_.end());
    ents_.push_back(id);
  } // insert


  //! Number of entities in this branch
  int size() const { return ents_.size(); }

  //! Remove a specific entity from the branch
  void remove(const size_t &id) {
    auto itr = find(ents_.begin(), ents_.end(), id);
    assert(itr != ents_.end());
    ents_.erase(itr);
  }


  auto begin() { return ents_.begin(); }
  auto end() { return ents_.end(); }
  auto clear() {
    ents_.clear();
    set_requested(false); 
    set_ghosts_local(false);
  }

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
  point_t bmin_;
  point_t bmax_;
  int owner_;
  point_t coordinates_;
  element_t mass_;
  std::vector<size_t> ents_;
  /** Bit representation of boolean value
   * | rqsted_ | ghs_loc | locality_ 2bits | leaf_ | children 8bits | 
   * Default rqsted = false, ghs_loc = true, loc = EMPTY, leaf = true, 8(0)
   */
  unsigned int type_ =  
      (uint(1)<<bits::LEAF) 
    | (uint(1)<<bits::GHS_LOC) 
    | (uint(EMPTY) << bits::LOC);
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
