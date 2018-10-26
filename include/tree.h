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
#include "utils.h"

#include "body.h"

using namespace flecsi;

namespace flecsi{
namespace execution{
void specialization_driver(int argc, char * argv[]);
void driver(int argc, char*argv[]); 
} // namespace execution
} // namespace flecsi 

struct body_holder_mpi_t{
  static const size_t dimension = gdimension;
  using element_t = type_t; 
  using point_t = flecsi::point__<element_t, dimension>;

  point_t position; 
  int owner; 
  double mass;
  flecsi::topology::entity_id_t id;
  double h; 
};

class tree_policy{
public:
  using tree_t = flecsi::topology::tree_topology<tree_policy>;
  using branch_int_t = uint64_t;
  static const size_t dimension = gdimension;
  using element_t = type_t; 
  using point_t = flecsi::point__<element_t, dimension>;
  using space_vector_t = flecsi::space_vector<element_t,dimension>;
  using geometry_t = flecsi::topology::tree_geometry<element_t, gdimension>;
  using id_t = flecsi::topology::entity_id_t;
  using key_t = flecsi::topology::morton_id<branch_int_t,dimension>;


  /** 
   * @brief BODY_HOLDER, entity in the tree. Light representation of a body
   */ 
  class body_holder : 
    public flecsi::topology::tree_entity<branch_int_t,dimension>{
  
  public: 

    body_holder(
        point_t coordinates,
        body * bodyptr,
        size_t owner,
        element_t mass,
	      id_t id,
        element_t h
        )
      : coordinates_(coordinates),mass_(mass),h_(h),bodyptr_(bodyptr)
    {
      locality_ = bodyptr_==nullptr?NONLOCAL:EXCL;
      global_id_ = id; 
      owner_ = owner;
    };

    body_holder()
      :coordinates_(point_t{0,0,0}),
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
    
    point_t coordinates() const {return coordinates_;}
    body* getBody(){return bodyptr_;};
    element_t mass(){return mass_;};
    element_t h(){return h_;};
    
    void setBody(body * bodyptr){bodyptr_ = bodyptr;};
    void set_id(id_t& id){id_ = id;}; 
    void set_coordinates(point_t& coordinates){coordinates_=coordinates;};
    void set_h(element_t h){h_=h;};

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
    point_t coordinates_;
    element_t mass_;
    element_t h_; 
    body * bodyptr_; 
  };
    
  using entity_t = body_holder;
   
  /**
   *  @brief BRANCH structure, a Center of Mass in our case.
   */ 
  class branch : public flecsi::topology::tree_branch<branch_int_t,dimension,
  double>{
  public:
    branch(){}

    branch(const branch_id_t& id):tree_branch(id){}

    void insert(body_holder* ent){
      ents_.push_back(ent); 
      if(ents_.size() > num_children){
        refine();
      }
    } // insert

    // Copy constructor 
    branch(const branch& b){
      this->ents_ = b.ents_;
      this->bmax_ = b.bmax_; 
      this->bmin_ = b.bmin_;
      this->mass_ = b.mass_; 
      this->radius_ = b.radius_;
      this->coordinates_ = b.coordinates_; 
      this->id_ = b.id_; 
      this->sub_entities_ = b.sub_entities_; 
      this->leaf_ = b.leaf_;  
    }

    ~branch(){
      ents_.clear();
    }
    
    auto begin(){
      return ents_.begin();
    }

    auto end(){
      return ents_.end();
    }

    auto clear(){
      ents_.clear(); 
    }

    void remove(body_holder* ent){
      auto itr = find(ents_.begin(), ents_.end(), ent);
      ents_.erase(itr);  
      if(ents_.empty()){
        coarsen();
      } 
    }

    point_t 
    coordinates(
        const std::array<flecsi::point__<element_t, dimension>,2>& range) 
    const{
      point_t p;
      branch_id_t bid = id(); 
      bid.coordinates(range,p);
      return p;
    }

    point_t coordinates(){return coordinates_;};
    element_t mass(){return mass_;};
    element_t radius(){return radius_;};
    point_t bmin(){return bmin_;};
    point_t bmax(){return bmax_;};
    void set_coordinates(point_t& coordinates){coordinates_=coordinates;};
    void set_mass(element_t mass){mass_ = mass;};
    void set_radius(element_t radius){radius_ = radius;};
    void set_bmax(point_t bmax){bmax_ = bmax;};
    void set_bmin(point_t bmin){bmin_ = bmin;};

   private:
    point_t coordinates_; 
    double mass_; 
    double radius_;
    std::vector<body_holder*> ents_; 
    point_t bmax_;
    point_t bmin_;
  }; // class branch 

  bool should_coarsen(branch* parent){
    return true;
  }

  using branch_t = branch;

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
    element_t mass = element_t(0); 
    point_t bmax{}; 
    point_t bmin{}; 
    element_t radius; 
    point_t coordinates = point_t{};
    uint64_t nchildren = 0; 
    for(size_t d = 0 ; d < dimension ; ++d){
      bmax[d] = -DBL_MAX; 
      bmin[d] = DBL_MAX; 
    }
    if(b->is_leaf()){
      for(auto child: *b)
      {
        if(local_only && !child->is_local()){
          continue; 
        }
        nchildren++; 
        element_t childmass = child->mass();
        for(size_t d = 0 ; d < dimension ; ++d)
        {
          bmax[d] = std::max(bmax[d],child->coordinates()[d]+epsilon+
              child->h());
          bmin[d] = std::min(bmin[d],child->coordinates()[d]-epsilon-
              child->h());
        }
        coordinates += childmass * child->coordinates(); 
        mass += childmass;  
      }
      if(mass > element_t(0))
        coordinates /= mass; 
    }else{
      for(int i = 0 ; i < (1<<dimension); ++i)
      {
        auto branch = tree->child(b,i);
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
    }
    b->set_sub_entities(nchildren);
    b->set_coordinates(coordinates);
    b->set_mass(mass);
    b->set_bmin(bmin);
    b->set_bmax(bmax);
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
    res = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[3]);
  return res;
}

#endif // tree_h
