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
#include <mpi.h>

#include "flecsi/geometry/point.h"

#include "tree_entity_id.h"
#include "key_id.h"
#include "tree_branch.h"
#include "tree_geometry.h"
#include "entity.h"

namespace flecsi {
namespace topology {
/*!
  Tree entity base class.
 */
template<
  typename TT,
  typename T,
  size_t D,
  typename ENT
>
class tree_entity{
public:

  using type_t = TT;
  using id_t = entity_id_t;
  using key_id_t = key_id__<T,D>;
  using point_t = point__<TT,D>;
  using range_t = std::array<point_t,2>;
  using entity_t = ENT;

protected:
  enum e_locality_ {LOCAL=0,NONLOCAL=1,SHARED=2,EXCL=3,GHOST=4};

public:

  tree_entity(
    const key_id_t& key,
    const point_t& coordinates,
    entity_t* entity_ptr,
    const int64_t owner,
    const type_t& mass,
    const entity_id_t& id,
    const type_t& radius
  ):key_id_(key),coordinates_(coordinates),entity_ptr_(entity_ptr),
  owner_(owner),mass_(mass),id_(id),radius_(radius)
  {
    locality_ = entity_ptr_==nullptr?NONLOCAL:EXCL;
    global_id_ = id;
  };

  tree_entity()
  : key_id_(key_id_t::null()),
  locality_(NONLOCAL)
  {
    owner_ = -1;
    global_id_ = {};
  }

  tree_entity(
    const key_id_t& key,
    const point_t& coordinates
  )
  : key_id_(key),
  coordinates_(coordinates),
  locality_(NONLOCAL){};

  key_id_t
  get_entity_key() const
  {
    return key_id_;
  }

  key_id_t
  key() const
  {
    return key_id_;
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
    return key_id_ != key_id_t::null();
  }

  /*!
   * Return true if the entity is local in this process
   */
  bool
  is_local() const
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    return owner_ == rank;
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

  point_t
  coordinates() const
  {
    return coordinates_;
  }

  void
  set_coordinates(
    point_t& coordinates)
  {
    coordinates_=coordinates;
  };

  type_t mass(){return mass_;};
  void set_h(type_t h){radius_=h;};
  type_t h(){return radius_;};

  void set_shared(){locality_ = SHARED;};

  friend std::ostream& operator<<(std::ostream& os, const tree_entity& b){
    os << std::setprecision(10);
    os << "Tree entity. Pos: " <<b.coordinates_ << " Mass: "<< b.mass_ << " ";
    if(b.locality_ == LOCAL || b.locality_ == EXCL || b.locality_ == SHARED)
    {
      os<< "LOCAL";
    }else{
      os << "NONLOCAL";
    }
    os << " owner: " << b.owner_;
    os << " id: " << b.id_;
    return os;
  }

  entity_t* getBody()
  { return entity_ptr_; }

  void setBody(entity_t* entity_ptr){
    entity_ptr_ = entity_ptr;
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
  set_entity_key_(
    key_id_t bid
  )
  {
    key_id_ = bid;
  }

  type_t mass_;
  type_t radius_;
  point_t coordinates_;
  key_id_t key_id_;
  entity_id_t id_;
  entity_id_t global_id_;
  e_locality_ locality_;
  int64_t owner_;

  entity_t* entity_ptr_;
};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_entity_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
