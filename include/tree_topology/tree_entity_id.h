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

#ifndef flecsi_topology_tree_entity_id_h
#define flecsi_topology_tree_entity_id_h

/*!
  \file tree_topology.h
  \authors nickm@lanl.gov
  \date Initial file creation: Apr 5, 2016
 */

/*
  Tree topology is a statically configured N-dimensional hashed tree for
  representing localized entities, e.g. particles. It stores entities in a
  configurable branch type. Inserting entities into a branch can cause that
  branch to be refined or coarsened correspondingly. A client of tree topology
  defines a policy which defines its branch and entity types and other
  compile-time parameters. Specializations can define a policy and default
  branch types which can then be specialized in a simpler fashion
  (see the basic_tree specialization).
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

namespace flecsi {
namespace topology {

  /*!
  All tree entities have an associated entity id of this type which is needed
  to interface with the index space.
 */
class entity_id_t{
public:
  entity_id_t()
  {}

  entity_id_t(
    const entity_id_t& id
  )
  : id_(id.id_)
  {}

  entity_id_t(
    size_t id
  )
  : id_(id)
  {}

  operator size_t() const
  {
    return id_;
  }

  entity_id_t&
  operator=(
    const entity_id_t& id
  )
  {
    id_ = id.id_;
    return *this;
  }

  size_t 
  value_()
  {
    return id_;
  }

  size_t
  index_space_index() const
  {
    return id_;
  }

private:
  size_t id_;
};

} // namespace topology
} // namespace flecsi

#endif  // flecsi_topology_tree_entity_id_h
