/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
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
 * @file tree.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Description of the tree policy used in the tree_topology from FleCSI.
 */

#ifndef tree_h
#define tree_h

#include <vector>

//#warning "CHANGE TO FLECSI ONE"
#include "tree_topology/tree_topology.h"
#include "tree_topology/tree_entity_id.h"
#include "flecsi/geometry/point.h"
#include "flecsi/geometry/space_vector.h"
//#include "utils.h"

#include "body.h"

using namespace flecsi;

namespace flecsi{
namespace execution{
void specialization_driver(int argc, char * argv[]);
void driver(int argc, char*argv[]);
} // namespace execution
} // namespace flecsi

class tree_policy{
public:

  using tree_t = flecsi::topology::tree_topology<tree_policy>;
  using key_int_t = uint64_t;
  static const size_t dimension = gdimension;
  using element_t = type_t;

  using key_t = flecsi::topology::key_id__<key_int_t,gdimension>;

  using point_t = flecsi::point__<element_t, dimension>;
  using space_vector_t = flecsi::space_vector<element_t,dimension>;
  using geometry_t = flecsi::topology::tree_geometry<element_t, gdimension>;
  using id_t = flecsi::topology::entity_id_t;

  /**
   * @brief BODY_HOLDER, entity in the tree. Light representation of a body
   */
  class body_holder :
    public flecsi::topology::tree_entity<double,key_int_t,dimension>{

  public:

    body_holder(
        const key_t& key,
        const point_t coordinates,
        body * bodyptr,
        size_t owner,
        element_t mass,
	      id_t id,
        element_t h
        )
      : tree_entity(key,coordinates),
      mass_(mass),
      h_(h),
      bodyptr_(bodyptr)
    {
      locality_ = bodyptr_==nullptr?NONLOCAL:EXCL;
      global_id_ = id;
      owner_ = owner;
    };

    body_holder()
      :tree_entity(key_id_t{},point_t{0,0,0}),
       mass_(0.0),
       h_(0.0),
       bodyptr_(nullptr)
    {
      locality_ = NONLOCAL;
      owner_ = -1;
      global_id_ = {};
    };

    ~body_holder()
    {
      //bodyptr_ = nullptr;
    }

    body* getBody(){return bodyptr_;};
    element_t mass(){return mass_;};
    element_t h(){return h_;};

    void setBody(body * bodyptr){bodyptr_ = bodyptr;};
    void set_id(id_t& id){id_ = id;};
    void set_h(element_t h){h_=h;};
    void set_shared(){locality_ = SHARED;};

    friend std::ostream& operator<<(std::ostream& os, const body_holder& b){
      os << std::setprecision(10);
      os << "Holder. Pos: " <<b.coordinates_ << " Mass: "<< b.mass_ << " ";
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

  private:
    element_t mass_;
    element_t h_;
    body * bodyptr_;
  };

  using tree_entity_t = body_holder;
  using branch_t = flecsi::topology::tree_branch<key_int_t,dimension,double>;
  using entity_t = body;

}; // class tree_policy

using tree_topology_t = flecsi::topology::tree_topology<tree_policy>;
using tree_geometry_t = flecsi::topology::tree_geometry<type_t,gdimension>;
using body_holder = tree_topology_t::body_holder;
using point_t = tree_topology_t::point_t;
using branch_t = tree_topology_t::branch_t;
using branch_id_t = tree_topology_t::branch_id_t;
using space_vector_t = tree_topology_t::space_vector_t;
using entity_key_t = tree_topology_t::key_t;
using entity_id_t = flecsi::topology::entity_id_t;

using range_t = std::array<point_t,2>;

inline
bool
operator==(
    const point_t& p1,
    const point_t& p2)
{
  for(size_t i=0;i<gdimension;++i)
    if(p1[i]!=p2[i])
      return false;
  return true;
}

inline
bool
operator!=(
    const point_t& p1,
    const point_t& p2)
{
  for(size_t i=0;i<gdimension;++i)
    if(p1[i]!=p2[i])
      return true;
  return false;
}

inline
point_t
operator+(
    const point_t& p,
    const double& val)
{
  point_t pr = p;
  for(size_t i=0;i<gdimension;++i)
    pr[i]+=val;
  return pr;
}

inline
point_t
operator-(
    const point_t& p,
    const double& val)
{
  point_t pr = p;
  for(size_t i=0;i<gdimension;++i)
    pr[i]-=val;
  return pr;
}

inline
bool
operator<(
    const point_t& p,
    const point_t& q)
{
  for(size_t i=0;i<gdimension;++i)
    if(p[i]>q[i])
      return false;
  return true;
}

inline
bool
operator>(
    const point_t& p,
    const point_t& q)
{
  for(size_t i=0;i<gdimension;++i)
    if(p[i]<q[i])
      return false;
  return true;
}

inline
point_t
operator*(
    const point_t& p,
    const point_t& q)
{
  point_t r = p;
  for(size_t i=0;i<gdimension;++i)
    r[i] *= q[i];
  return r;
}

inline double norm_point( const point_t& p) {
  double res = 0;
  if constexpr (gdimension == 1)
    res = std::abs(p[0]);
  else if constexpr (gdimension == 2)
    res = sqrt(p[0]*p[0] + p[1]*p[1]);
  else
    res = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  return res;
}

#endif // tree_h
