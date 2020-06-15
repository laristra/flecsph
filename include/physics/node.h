/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2020 Triad National Security, LLC
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
 * @file cofm.h
 * @author Julien Loiseau
 * @date January 2020
 * @brief Representation of a center of mass
 */

#pragma once

#include "space_vector.h"
#include "tensor.h"
#include "tree_topology/tree_types.h"
#include "user.h"

namespace flecsi {
using sym_tensor_rank2 =
  flecsi::tensor_u<type_t, symmetry_type::symmetric, gdimension, gdimension>;
using sym_tensor_rank3 = flecsi::tensor_u<type_t,
  symmetry_type::symmetric,
  gdimension,
  gdimension,
  gdimension>;
using sym_tensor_rank4 = flecsi::tensor_u<type_t,
  symmetry_type::symmetric,
  gdimension,
  gdimension,
  gdimension,
  gdimension>;
} // namespace flecsi

template<class KEY, size_t FMM_ORDER>
class node_u : public flecsi::topology::cofm_u<gdimension, type_t, KEY>
{

  // unspecialized

}; // class node

//
// node_u: partial specialization for 1st-order Taylor expansion
//
template<class KEY>
class node_u<KEY, 1> : public flecsi::topology::cofm_u<gdimension, type_t, KEY>
{

  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::space_vector_u<element_t, dimension>;
  using key_t = KEY;

public:
  node_u() : flecsi::topology::cofm_u<gdimension, type_t, KEY>() {
    pc_ = 0;
    fc_ = 0;
    affected_ = false;
  }

  node_u(const key_t & key)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(key) {
    pc_ = 0;
    fc_ = 0;
    affected_ = false;
  }

  explicit node_u(const node_u & c)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(c) {
    pc_ = c.pc_;
    fc_ = c.fc_;
    affected_ = c.affected_;
  }

  const type_t & pc() const {
    return pc_;
  }
  const point_t & fc() const {
    return fc_;
  }

  type_t & pc() {
    return pc_;
  }
  point_t & fc() {
    return fc_;
  }

  void set_affected(const bool & affected) {
    affected_ = affected;
  }
  bool affected() const {
    return affected_;
  }

private:
  type_t pc_;
  point_t fc_;

  bool affected_;

}; // class node<KEY, 1>

//
// node_u: partial specialization for 2nd order Taylor expansion
//
template<class KEY>
class node_u<KEY, 2> : public flecsi::topology::cofm_u<gdimension, type_t, KEY>
{

  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::space_vector_u<element_t, dimension>;
  using key_t = KEY;

  using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
  using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
  using sym_tensor_rank4 = flecsi::sym_tensor_rank4;

public:
  node_u() : flecsi::topology::cofm_u<gdimension, type_t, KEY>() {
    Q_ = {0};
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    affected_ = false;
  }

  node_u(const key_t & key)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(key) {
    Q_ = {0};
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    affected_ = false;
  }

  explicit node_u(const node_u & c)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(c) {
    Q_ = c.quad();
    pc_ = c.pc_;
    fc_ = c.fc_;
    dfcdr_ = c.dfcdr_;
    affected_ = c.affected_;
  }

  const sym_tensor_rank2 & quad() const {
    return Q_;
  }

  sym_tensor_rank2 & quad() {
    return Q_;
  }

  const type_t & pc() const {
    return pc_;
  }
  const point_t & fc() const {
    return fc_;
  }
  const sym_tensor_rank2 & dfcdr() const {
    return dfcdr_;
  }

  type_t & pc() {
    return pc_;
  }
  point_t & fc() {
    return fc_;
  }
  sym_tensor_rank2 & dfcdr() {
    return dfcdr_;
  }

  void set_affected(const bool & affected) {
    affected_ = affected;
  }
  bool affected() const {
    return affected_;
  }

private:
  sym_tensor_rank2 Q_;

  type_t pc_;
  point_t fc_;
  sym_tensor_rank2 dfcdr_;

  bool affected_;

}; // class node<KEY, 2>

//
// node_u: partial specialization for 3rd order Taylor expansion
//
template<class KEY>
class node_u<KEY, 3> : public flecsi::topology::cofm_u<gdimension, type_t, KEY>
{

  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::space_vector_u<element_t, dimension>;
  using key_t = KEY;

  using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
  using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
  using sym_tensor_rank4 = flecsi::sym_tensor_rank4;

public:
  node_u() : flecsi::topology::cofm_u<gdimension, type_t, KEY>() {
    H_ = {0};
    Q_ = {0};
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    dfcdrdr_ = 0;
    affected_ = false;
  }

  node_u(const key_t & key)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(key) {
    H_ = {0};
    Q_ = {0};
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    dfcdrdr_ = 0;
    affected_ = false;
  }

  explicit node_u(const node_u & c)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(c) {
    H_ = c.octo();
    Q_ = c.quad();
    pc_ = c.pc_;
    fc_ = c.fc_;
    dfcdr_ = c.dfcdr_;
    dfcdrdr_ = c.dfcdrdr_;
    affected_ = c.affected_;
  }

  const sym_tensor_rank3 & octo() const {
    return H_;
  }
  const sym_tensor_rank2 & quad() const {
    return Q_;
  }

  sym_tensor_rank3 & octo() {
    return H_;
  }
  sym_tensor_rank2 & quad() {
    return Q_;
  }

  const type_t & pc() const {
    return pc_;
  }
  const point_t & fc() const {
    return fc_;
  }
  const sym_tensor_rank2 & dfcdr() const {
    return dfcdr_;
  }
  const sym_tensor_rank3 & dfcdrdr() const {
    return dfcdrdr_;
  }

  type_t & pc() {
    return pc_;
  }
  point_t & fc() {
    return fc_;
  }
  sym_tensor_rank2 & dfcdr() {
    return dfcdr_;
  }
  sym_tensor_rank3 & dfcdrdr() {
    return dfcdrdr_;
  }

  void set_affected(const bool & affected) {
    affected_ = affected;
  }
  bool affected() const {
    return affected_;
  }

private:
  sym_tensor_rank3 H_;
  sym_tensor_rank2 Q_;

  type_t pc_;
  point_t fc_;
  sym_tensor_rank2 dfcdr_;
  sym_tensor_rank3 dfcdrdr_;

  bool affected_;

}; // class node<KEY, 3>

//
// node_u: partial specialization for 4rd order Taylor expansion
//
template<class KEY>
class node_u<KEY, 4> : public flecsi::topology::cofm_u<gdimension, type_t, KEY>
{

  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::space_vector_u<element_t, dimension>;
  using key_t = KEY;

  using sym_tensor_rank2 = flecsi::sym_tensor_rank2;
  using sym_tensor_rank3 = flecsi::sym_tensor_rank3;
  using sym_tensor_rank4 = flecsi::sym_tensor_rank4;

public:
  node_u() : flecsi::topology::cofm_u<gdimension, type_t, KEY>() {
    X_ = {0};
    H_ = {0};
    Q_ = {0};
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    dfcdrdr_ = 0;
    dfcdrdrdr_ = 0;
    affected_ = false;
  }

  node_u(const key_t & key)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(key) {
    X_ = {0};
    H_ = {0};
    Q_ = {0};
    pc_ = 0;
    fc_ = 0;
    dfcdr_ = 0;
    dfcdrdr_ = 0;
    dfcdrdrdr_ = 0;
    affected_ = false;
  }

  explicit node_u(const node_u & c)
    : flecsi::topology::cofm_u<gdimension, type_t, KEY>(c) {
    X_ = c.hexa();
    H_ = c.octo();
    Q_ = c.quad();
    pc_ = c.pc_;
    fc_ = c.fc_;
    dfcdr_ = c.dfcdr_;
    dfcdrdr_ = c.dfcdrdr_;
    dfcdrdrdr_ = c.dfcdrdrdr_;
    affected_ = c.affected_;
  }

  const sym_tensor_rank4 & hexa() const {
    return X_;
  }
  const sym_tensor_rank3 & octo() const {
    return H_;
  }
  const sym_tensor_rank2 & quad() const {
    return Q_;
  }

  sym_tensor_rank4 & hexa() {
    return X_;
  }
  sym_tensor_rank3 & octo() {
    return H_;
  }
  sym_tensor_rank2 & quad() {
    return Q_;
  }

  const type_t & pc() const {
    return pc_;
  }
  const point_t & fc() const {
    return fc_;
  }
  const sym_tensor_rank2 & dfcdr() const {
    return dfcdr_;
  }
  const sym_tensor_rank3 & dfcdrdr() const {
    return dfcdrdr_;
  }
  const sym_tensor_rank4 & dfcdrdrdr() const {
    return dfcdrdrdr_;
  }

  type_t & pc() {
    return pc_;
  }
  point_t & fc() {
    return fc_;
  }
  sym_tensor_rank2 & dfcdr() {
    return dfcdr_;
  }
  sym_tensor_rank3 & dfcdrdr() {
    return dfcdrdr_;
  }
  sym_tensor_rank4 & dfcdrdrdr() {
    return dfcdrdrdr_;
  }

  void set_affected(const bool & affected) {
    affected_ = affected;
  }
  bool affected() const {
    return affected_;
  }

private:
  sym_tensor_rank4 X_;
  sym_tensor_rank3 H_;
  sym_tensor_rank2 Q_;

  type_t pc_;
  point_t fc_;
  sym_tensor_rank2 dfcdr_;
  sym_tensor_rank3 dfcdrdr_;
  sym_tensor_rank4 dfcdrdrdr_;

  bool affected_;

}; // class node<KEY, 4>
