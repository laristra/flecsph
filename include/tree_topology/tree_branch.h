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

#include "morton_id.h"

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
  using branch_id_t = morton_id<T, D>;
  using id_t = branch_id_t;
  static constexpr size_t num_children = 1 << dimension;
  using point_t = point__<E,D>;
  using element_t = E;

  tree_branch()
  : action_(action::none)
  {}

  tree_branch(const branch_id_t& id)
  : action_(action::none),
  id_(id)
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
    return leaf_;
  }

  void 
  set_leaf(
      bool leaf)
  {
    leaf_ = leaf; 
  }

  bool
  is_valid() const
  {
    return true;
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
 
  action action_;
  branch_id_t id_;
  uint64_t sub_entities_ = 0; // Subentities in this subtree
  bool leaf_ = true; 

};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_branch_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
