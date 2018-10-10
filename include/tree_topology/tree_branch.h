/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_topology_tree_branch_h
#define flecsi_topology_tree_branch_h

/*!
  \file tree_branch.h
  \authors jloiseau@lanl.gov
  \date Initial file creation: Oct 9 2018
 */

#include <map>
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <functional>
#include <mutex>
#include <stack>
#include <math.h>
#include <float.h>

#include "flecsi/geometry/point.h"
#include "flecsi/concurrency/thread_pool.h"
#include "flecsi/data/storage.h"
#include "flecsi/data/data_client.h"
#include "flecsi/topology/index_space.h"

#include "morton_branch_id.h"

/*
#define np(X)                                                            \
 std::cout << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ \
           << ": " << #X << " = " << (X) << std::endl

#define hp(X)                                                            \
 std::cout << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ \
           << ": " << #X << " = " << std::hex << (X) << std::endl
*/

namespace flecsi {
namespace topology {


/*! When an entoty is added or removed from a branch, the user-level tree 
 * may trigger one of these actions. 
 */
enum class action: uint8_t{
  none=0b00, 
  refine = 0b01, 
  coarsen = 0b10
};


/*!
  Tree branch base class.
 */
template<
  typename T,
  size_t D,
  typename E
>
class tree_branch
{
public:
  using branch_int_t = T;

  static const size_t dimension = D;

  using branch_id_t = branch_id<T, D>;

  using id_t = branch_id_t;

  static constexpr size_t num_children = branch_int_t(1) << dimension;

  using point_t = point__<E,D>;
  using element_t = E;

  tree_branch()
  : action_(action::none),
  coordinates_(point_t{})
  {}

  tree_branch(const branch_id_t& id)
  : action_(action::none),
  id_(id),
  coordinates_(point_t{})
  {}

  branch_id_t
  id() const
  {
    return id_;
  }

  /*!
    Called to trigger a refinement at this branch.
   */
  void refine()
  {
    action_ = action::refine;
  }

  /*!
    Called to trigger a coarsening at this branch.
   */
  void coarsen()
  {
    action_ = action::coarsen;
  }

  /*!
    Clear refine/coarsen actions.
   */
  void reset()
  {
    action_ = action::none;
  }

  bool
  is_leaf() const
  {
    return nchildren_!=0;
  }

  bool
  is_valid() const
  {
    return true;
  }

  void 
  set_mass(element_t mass)
  {
    mass_ = mass;
  }
 
  void 
  set_coordinates(point_t coordinates)
  {
    coordinates_ = coordinates;
  }

  element_t 
  mass() const
  {
    return mass_;
  }

  point_t 
  bmin() const 
  {
    return bmin_;
  }

  point_t
  bmax() const 
  {
    return bmax_;
  }

  void 
  set_bmin(
    point_t bmin)
  {
    bmin_ = bmin;
  }
 
  void 
  set_bmax(
    point_t bmax)
  {
    bmax_ = bmax;
  }

  point_t 
  get_coordinates() const
  {
    return coordinates_;
  }

  uint64_t 
  sub_entities() const
  {
    return sub_entities_;
  }

  void 
  set_sub_entities(
      uint64_t sub_entities)
  {
    sub_entities_ = sub_entities;
  }


protected:
  template<class P>
  friend class tree_topology;

  void
  set_id_(
    branch_id_t id
  )
  {
    id_ = id;
  }

  action
  requested_action_()
  {
    return action_;
  }

  /*template<
    class B
  >
  B*
  child_(
    size_t ci
  ) const
  {
    assert(ci < num_children);
    return static_cast<B*>(children_) + ci;
  }

  template<
    class B
    >
  bool
  into_branch_()
  {
    if(children_)
    {
      return false;
    }
    auto c = new B[num_children];
    for(branch_int_t bi = 0; bi < num_children; ++bi)
    {
      B& ci = c[bi];
      ci.id_ = id_;
      ci.id_.push(bi);
    }
    children_ = c;
    return true;
  }

  template<
    class B
  >
  void
  into_leaf_()
  {
    if(children_)
    {
      delete[] static_cast<B*>(children_);
      children_ = nullptr;
    }
  }

  template<
    class B
  >
  void
  dealloc_()
  {
    if(children_)
    {
      for(size_t i = 0; i < num_children; ++i)
      {
        static_cast<B*>(children_)[i].template dealloc_<B>();
      }

      delete[] static_cast<B*>(children_);
      children_ = nullptr;
    }
  }*/

  action action_;

  branch_id_t id_;

  point_t coordinates_; 
  //element_t radius_; 
  element_t mass_;

  point_t bmin_; 
  point_t bmax_; 

  uint64_t nchildren_; // Number of bodies attached
  uint64_t sub_entities_; // Subentities in the leaves the subtree
};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_branch_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
