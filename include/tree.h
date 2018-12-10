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
      bodyptr_ = nullptr;
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

template <
  typename T,
  size_t D
>
struct traversal
{
  using element_t = T;
  const static size_t dimension = D;


 // Functions for the tree traversal
  static void update_COM(
      tree_topology_t* tree,
      branch_t * b,
      element_t epsilon = element_t(0),
      bool local_only = false)
  {
    bool local_branch = false;
    element_t mass = element_t(0);
    point_t bmax{};
    point_t bmin{};
    element_t radius = 0.;
    point_t coordinates = point_t{};
    uint64_t nchildren = 0;
    for(size_t d = 0 ; d < dimension ; ++d){
      bmax[d] = -DBL_MAX;
      bmin[d] = DBL_MAX;
    }
    if(b->is_leaf()){
      // For local branches, compute the radius
      if(b->is_local()){
        for(auto child: *b)
        {
          auto ent = tree->get(child);
          if(ent->is_local()){
            local_branch = true;
          }
          if(local_only && !ent->is_local()){
            continue;
          }
          nchildren++;
          element_t childmass = ent->mass();
          for(size_t d = 0 ; d < dimension ; ++d)
          {
            bmax[d] = std::max(bmax[d],ent->coordinates()[d]+epsilon+ent->h());
            bmin[d] = std::min(bmin[d],ent->coordinates()[d]-epsilon-ent->h());
          }
          coordinates += childmass * ent->coordinates();
          mass += childmass;
        }
        if(mass > element_t(0))
          coordinates /= mass;
        // Compute the radius
        for(auto child: *b)
        {
          auto ent = tree->get(child);
          radius = std::max(
              radius,
              distance(ent->coordinates(),coordinates)  + epsilon + ent->h());
        }
      }else{
        // For non local particles use existing value from remote 
        coordinates = b->coordinates();
        bmin = b->bmin();
        bmax = b->bmax();
        mass = b->mass();
        nchildren = 1;
        local_branch = false;
      }
    }else{
      for(int i = 0 ; i < (1<<dimension); ++i)
      {
        auto branch = tree->child(b,i);
        if(branch->is_local()){
          local_branch=true;
        }
        nchildren+=branch->sub_entities();
        mass += branch->mass();
        if(branch->mass() > 0)
        {
          for(size_t d = 0 ; d < dimension ; ++d)
          {
            bmax[d] = std::max(bmax[d],branch->bmax()[d]);
            bmin[d] = std::min(bmin[d],branch->bmin()[d]);
          }
        }
        coordinates += branch->mass()*branch->coordinates();
      }
      if(mass > element_t(0))
        coordinates /= mass;
      // Compute the radius
      for(int i = 0 ; i < (1<<dimension); ++i)
      {
        auto branch = tree->child(b,i);
        radius = std::max(
            radius,
            distance(coordinates,branch->coordinates()) + branch->radius());
      }
    }
    b->set_radius(radius);
    b->set_sub_entities(nchildren);
    b->set_coordinates(coordinates);
    b->set_mass(mass);
    b->set_bmin(bmin);
    b->set_bmax(bmax);
    local_branch?
      b->set_locality(branch_t::LOCAL):
      b->set_locality(branch_t::NONLOCAL);
    if(!b->is_local()) tree->nonlocal_branches_add();
  }


};

using traversal_t = traversal<double,gdimension>;

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
