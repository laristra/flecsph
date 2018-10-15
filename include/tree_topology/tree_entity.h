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
  enum locality {LOCAL=0,NONLOCAL=1,SHARED=2,EXCL=3,GHOST=4}; 
  
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
  setLocality(locality loc)
  {
    locality_ = loc;
  }
  
  locality 
  getLocality()
  {
    return locality_;
  };

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
  set_branch_id_(
    branch_id_t bid
  )
  {
    branch_id_ = bid;
  }

  branch_id_t branch_id_;
  entity_id_t id_;

  locality locality_;
};
 
} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_entity_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
