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

#ifndef flecsi_topology_tree_entity_h
#define flecsi_topology_tree_entity_h

/*!
  \file tree_entity.h
  \authors nickm@lanl.gov
  \date Initial file creation: Apr 5, 2016
 */


#include <map>
#include <unordered_map>
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

#include "tree_entity_id.h"
#include "morton_id.h"
#include "tree_branch.h"
#include "tree_geometry.h"

namespace flecsi {
namespace topology {

  
/*!
  Tree entity base class.
 */
template<
  typename T,
  size_t D
>
class tree_entity{
public:

  using id_t = entity_id_t;
  using branch_id_t = morton_id<T, D>;
  
protected:
  enum e_locality_ {LOCAL=0,NONLOCAL=1,SHARED=2,EXCL=3,GHOST=4}; 
  
public:

  tree_entity()
  : branch_id_(branch_id_t::null()),
  locality_(NONLOCAL)
  {}

  branch_id_t
  get_branch_id() const
  {
    return branch_id_;
  }

  entity_id_t
  id() const
  {
    return id_;
  }

  entity_id_t 
  global_id() const 
  {
    return global_id_;
  }

  void 
  set_global_id(entity_id_t id)
  {
    global_id_ = id;
  }

  entity_id_t
  index_space_id() const
  {
    return id_;
  }

  /*!
    Return whether the entity is current inserted in a tree.
   */
  bool
  is_valid() const
  {
    return branch_id_ != branch_id_t::null();
  }

  /*!
   * Return true if the entity is local in this process
   */
  bool
  is_local() const 
  {
    return (locality_ == LOCAL || locality_ == EXCL || locality_ == SHARED); 
  }

  void 
  set_locality(e_locality_ loc)
  {
    locality_ = loc;
  }
  
  e_locality_ 
  locality()
  {
    return locality_;
  }

  void 
  set_owner(int64_t owner)
  {
    owner_ = owner;
  }

  int64_t 
  owner()
  {
    return owner_; 
  }

protected:
  template<class P>
  friend class tree_topology;

  void
  set_id_(
    entity_id_t id
  )
  {
    id_ = id;
  }

  void
  set_global_id_(
    entity_id_t id
  )
  {
    global_id_ = id;
  }

  void
  set_branch_id_(
    branch_id_t bid
  )
  {
    branch_id_ = bid;
  }

  branch_id_t branch_id_;
  entity_id_t id_;
  entity_id_t global_id_;
  e_locality_ locality_;
  int64_t owner_;
};
 
} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_entity_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
