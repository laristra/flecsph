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

#ifndef flecsi_topology_tree_topology_h
#define flecsi_topology_tree_topology_h

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

#include "morton_id.h"
#include "tree_branch.h"
#include "tree_entity.h"
#include "tree_geometry.h"

namespace flecsi {
namespace topology {


// Hasher for the branch id used in the unordered_map data structure
template<
  typename T,
  size_t D,
  class IDTYPE=morton_id<T,D>
>
struct branch_id_hasher__{
  size_t
  operator()(
    const IDTYPE& k
  ) const
  {
    return std::hash<T>()(k.value_());
  }
};

/*!
  The tree topology is parameterized on a policy P which defines its branch and
  entity types.
 */
template<
  class P
>
class tree_topology : public P, public data::data_client_t
{
public:
  using Policy = P;

  static const size_t dimension = Policy::dimension;
  using element_t = typename Policy::element_t;
  using point_t = point__<element_t, dimension>;
  using range_t = std::pair<element_t, element_t>;
  using branch_int_t = typename Policy::branch_int_t;
  using branch_id_t = typename Policy::key_t; 
  using branch_id_vector_t = std::vector<branch_id_t>;

  using branch_t = typename Policy::branch_t;
  using branch_vector_t = std::vector<branch_t*>;

  using entity_t = typename Policy::entity_t;
  using entity_vector_t = std::vector<entity_t*>;
  using apply_function = std::function<void(branch_t&)>;
  using entity_id_vector_t = std::vector<entity_id_t>;
  using geometry_t = tree_geometry<element_t, dimension>;
  using entity_space_t = std::vector<entity_t>;
  using entity_space_ptr_t = std::vector<entity_t*>;

  struct filter_valid{
    bool operator()(entity_t* ent) const{
      return ent->is_valid();
    }
  };

  /*!
    Constuct a tree topology with unit coordinates, i.e. each coordinate
    dimension is in range [0, 1].
   */
  tree_topology()
  {
    branch_map_.emplace(branch_id_t::root(), branch_id_t::root());
    root_ = branch_map_.find(branch_id_t::root()); 
    assert(root_ != branch_map_.end());

    max_depth_ = 0;
    max_scale_ = element_t(1);

    for(size_t d = 0; d < dimension; ++d)
    {
      range_[0][d] = element_t(0);
      range_[1][d] = element_t(1);
      scale_[d] = element_t(1);
    }
  }

  /*!
    Construct a tree topology with specified ranges [end, start] for each
    dimension.
   */
  tree_topology(
    const point__<element_t, dimension>& start,
    const point__<element_t, dimension>& end
  )
  {
    branch_map_.emplace(branch_id_t::root(),branch_id_t::root());
    root_ = branch_map_.find(branch_id_t::root()); 
    assert(root_ != branch_map_.end());
   
    max_depth_ = 0;
    max_scale_ = end[0] - start[0]; 

    for(size_t d = 0; d < dimension; ++d)
    {
      scale_[d] = end[d] - start[d];
      max_scale_ = std::max(max_scale_, scale_[d]);
      range_[0][d] = start[d];
      range_[1][d] = end[d];
    }
  }

  /** 
   * @brief Destroy the tree: empty the hash-table and destroy the entities 
   * lists 
   */
  ~tree_topology()
  {
    branch_map_.clear();  
    entities_.clear(); 
    //entities_ghosts_.clear(); 
    ghosts_id_.clear(); 
  } 

  /** 
   * @brief Get the range of the current 
   */
  const std::array<point__<element_t,dimension>,2>& range()
  {
    return range_; 
  }

  /**
   * @brief Get the ci-th child of the given branch.
   */
  branch_t*
  child(
    branch_t* b,
    size_t ci
  )
  {
    // Use the hash table 
    branch_id_t bid = b->id(); // Branch id 
    bid.push(ci); // Add child number 
    auto child = branch_map_.find(bid); // Search for the child
    // If it does not exists, return nullptr 
    if(child == branch_map_.end())
    {
      return nullptr; 
    } 
    return &child->second;
  }

  /*!
    Return an index space containing all entities (including those removed).
   */
  std::vector<entity_t>&
  all_entities() const
  {
    return entities_;
  }

  /*!
    Return an index space containing all non-removed entities.
   */
  std::vector<entity_t>&
  entities()
  {
    return entities_;
  }

  //std::vector<entity_id_t>& 
  //entities_ghosts()
  //{
  //  return entities_ghosts_; 
  //}

  entity_t* 
  get_ghost(
      const size_t& global_id)
  {
    auto local_id_itr = ghosts_id_.find(global_id);
    assert(local_id_itr != ghosts_id_.end());
    auto local_id = local_id_itr->second; 
    return &(entities_[local_id]); 
  }

  /*!
    Insert an entity into the lowest possible branch division.
   */
  void
  insert(
    const entity_id_t& id
  )
  {
    insert(id, max_depth_);
  }

  /**
   * @brief Generic postorder traversal  
   * First version using two stacks.
   * \TODO switch to one stack  
   */
  template<
    typename F, 
    typename... ARGS
  >
  void 
  post_order_traversal(
      branch_t * start,
      F&& function,
      ARGS&&... args)
  {
    nonlocal_branches_ = 0L; 
    std::stack<branch_t*> stk1; 
    std::stack<branch_t*> stk2;
    stk1.push(start);
    while(!stk1.empty())
    {
      branch_t * cur = stk1.top(); 
      stk1.pop();
      stk2.push(cur);
      // Push children to stk1 
      if(!cur->is_leaf()){
        for(int i = 0 ; i < (1<<dimension) ; ++i)
        {
          branch_t * next = child(cur,i);
          stk1.push(next);
        }
      }
    }
    // Unstack stack 2 
    while(!stk2.empty())
    {
      branch_t * cur = stk2.top(); 
      stk2.pop();
      function(this,cur,std::forward<ARGS>(args)...);
    }
  }
 
  void
  get_all_branches(
    branch_t * start,
    std::vector<branch_t*>& search_list)
  {
    std::stack<branch_t*> stk;
    stk.push(start);
    search_list.push_back(start);
    
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        //search_list.push_back(c);
      }else{
        for(int i=0; i<(1<<dimension);++i){
          branch_t * next = child(c,i);
          if(next->sub_entities() > 0){
            search_list.push_back(next);
            stk.push(next);
          }
        }
      } 
    }
  }

  /**
   * @brief Return a vector with all the local sub entities
   * 
   * @param start The branch in which the search occur 
   * @param search_list The found entities vector
   */
  void
  get_sub_entities_local(
    branch_t * start,
    std::vector<entity_t*>& search_list)
  {
    std::stack<branch_t*> stk;
    stk.push(start);
    
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        for(auto bh: *c){
          if(bh->is_local()){
            search_list.push_back(bh);
          }
        }
      }else{
        for(int i=0; i<(1<<dimension);++i){
          branch_t * next = child(c,i);
          if(next->sub_entities() > 0){
            stk.push(next);
          }
        }
      } 
    }
  }

  void 
  find_sub_cells(
      branch_t * b,
      uint64_t criterion,
      std::vector<branch_t*>& search_list)
  {
    std::stack<branch_t*> stk;
    stk.push(b);
    
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf() && c->sub_entities() > 0){
        search_list.push_back(c);
      }else{
        if(c->sub_entities() <= criterion && c->sub_entities() > 0){
          search_list.push_back(c); 
        }else{
          for(int i=0; i<(1<<dimension);++i){
            branch_t * next = child(c,i);
            if(next->sub_entities() > 0){
              stk.push(next);
            }
          }
        }
      } 
    }
  }

  /**
   * @brief Find all the center of mass of the tree up to the 
   * maximum mass criterion.  
   * 
   * @param b The starting branch for the search, usually root
   * @param mass_criterion The maximum mass for the COMs
   * @param search_list The extracted COMs
   */
  void 
  find_sub_cells_mass(
      branch_t * b,
      double mass_criterion,
      std::vector<branch_t*>& search_list)
  {
    std::stack<branch_t*> stk;
    stk.push(b);
    
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf() && c->sub_entities() > 0){
        search_list.push_back(c);
      }else{
        if(c->mass() <= mass_criterion && c->sub_entities() > 0){
          search_list.push_back(c); 
        }else{
          for(int i=0; i<(1<<dimension);++i){
            branch_t * next = child(c,i);
            if(next->sub_entities() > 0){
              stk.push(next);
            }
          }
        }
      } 
    }
  }
  
  template<
    typename EF,
    typename... ARGS
  >
  void 
  apply_sub_cells(
      branch_t * b,
      element_t MAC,
      int64_t ncritical,
      bool do_square,
      EF&& ef,
      ARGS&&... args)
  {
    std::stack<branch_t*> stk;
    stk.push(b);

    std::vector<branch_t*> work_branch; 
    
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf() && c->sub_entities() > 0){
        work_branch.push_back(c);
      }else{
        if((int64_t)c->sub_entities() < ncritical && c->sub_entities() > 0){
          work_branch.push_back(c);
        }else{
          for(int i=0; i<(1<<dimension);++i){
            branch_t * next = child(c,i);
            if(next->sub_entities() > 0 && next->is_local()){
              stk.push(next);
            }
          }
        }
      } 
    }

    // Start the threads on the branches 
    #pragma omp parallel for 
    for(size_t i = 0; i < work_branch.size(); ++i){
      std::vector<branch_t*> inter_list; 
      sub_cells_inter(work_branch[i],MAC,inter_list);
      force_calc(work_branch[i],inter_list,do_square,
               ef,std::forward<ARGS>(args)...);
    }

  }

  void
  sub_cells_inter(
    branch_t* b,
    element_t MAC,
    std::vector<branch_t*>& inter_list)
  {
    std::stack<branch_t*> stk; 
    stk.push(root());
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        inter_list.push_back(c);
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(c,i);
          if(branch->sub_entities() > 0 && geometry_t::intersects_box_box(
            b->bmin(),
            b->bmax(),
            branch->bmin(),
            branch->bmax()))
          {
            stk.push(branch);
          }
        }
      }
    }
  }

  template<
    typename EF,
    typename... ARGS
  >
  void
  force_calc(
      branch_t* b, 
      std::vector<branch_t*>& inter_list,
      bool do_square, 
      EF&& ef,
      ARGS&&... args)
  {
    std::stack<branch_t*> stk;
    stk.push(b);
    while(!stk.empty())
    {
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        for(auto id: *c){
          auto child = this->get(id); 
          if(child->is_local()){
            if(do_square)
              apply_sub_entity_sq(child,inter_list,ef,
                std::forward<ARGS>(args)...);
            else 
              apply_sub_entity(child,inter_list,ef,
                std::forward<ARGS>(args)...);
          }
        }
      }else{
        for(int i=0; i<(1<<dimension); ++i){
          branch_t * next = child(c,i);
          if(next->mass() > 0.){
            stk.push(next);
          }
        }
      }
    }
  }

  template<
    typename EF,
    typename... ARGS
  >
  void 
  apply_sub_entity(
      entity_t* ent, 
      std::vector<branch_t*>& inter_list,
      EF&& ef,
      ARGS&&... args)
  {
    std::vector<entity_t*> neighbors;
    for(auto b: inter_list){
      for(auto id: *b){
        auto nb = this->get(id);
        if(geometry_t::within(
              ent->coordinates(),
              nb->coordinates(),
              ent->h()))
        {
          neighbors.push_back(nb);
        }
      }
    }
    ef(ent,neighbors,std::forward<ARGS>(args)...);
  }

  template<
    typename EF,
    typename... ARGS
  >
  void 
  apply_sub_entity_sq(
      entity_t* ent, 
      std::vector<branch_t*>& inter_list,
      EF&& ef,
      ARGS&&... args)
  {
    std::vector<entity_t*> neighbors;
    for(auto b: inter_list){
      for(auto id: *b){
        auto nb = this->get(id); 
        if(geometry_t::within_square(
              ent->coordinates(),
              nb->coordinates(),
              ent->h(),nb->h()))
        {
          neighbors.push_back(nb);
        }
      }
    }
    ef(ent,neighbors,std::forward<ARGS>(args)...);
  }


  /*!
    Return an index space containing all entities within the specified
    spheroid.
   */
  template<
    typename EF>
  entity_space_ptr_t
  find_in_radius(
    const point_t& center,
    element_t radius,
    EF&& ef
  )
  {
    entity_space_ptr_t ents;

    // ITERATIVE VERSION
    std::stack<branch_t*> stk;
    stk.push(root());

    while(!stk.empty()){
      branch_t* b = stk.top();
      stk.pop();
      if(b->is_leaf()){
        for(auto id: *b){
          auto child = &(entities_[id]);
          // Check if in radius 
          if(ef(center,child->coordinates(),radius,child->h())){
            ents.push_back(child);
          }
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(b,i);
          if(geometry_t::intersects_sphere_box(
                branch->bmin(),
                branch->bmax(),
                center,
                radius
                ))
          {
            stk.push(branch);
          }
        }
      }
    }
    return ents;
  }


/*!
    Return an index space containing all entities within the specified
    Box
   */
  template<
    typename EF> 
  entity_space_ptr_t
  find_in_box(
    const point_t& min,
    const point_t& max,
    EF&& ef
  )
  {
    entity_space_ptr_t ents;

    // ITERATIVE VERSION
    std::stack<branch_t*> stk;
    stk.push(root());

    while(!stk.empty()){
      branch_t* b = stk.top();
      stk.pop();
      if(b->is_leaf()){
        for(auto id: *b){
          auto child = &(entities_[id]);
          // Check if in box 
          if(ef(min,max,child->coordinates(),child->h())){
            ents.push_back(child);
          }
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(b,i);
          if(geometry_t::intersects_box_box(min,max,branch->bmin(),
                branch->bmax()))
          {
            stk.push(branch);
          }
        }
      }
    }
    return ents;
  }


    /*!
      Construct a new entity. The entity's constructor should not be called
      directly.
     */
    template<
      class... Args
    >
    entity_id_t
    make_entity(
      Args&&... args
    )
    {
      entities_.emplace_back(std::forward<Args>(args)...);
      auto ent = &(entities_.back());
      // Size -1 to start at 0 
      entity_id_t id = entities_.size()-1;
      ent->set_id_(id);
      if(!ent->is_local()){
        //entities_ghosts_.push_back(id);
        // Map of local-global id 
        ghosts_id_.insert(std::make_pair(ent->global_id(),id));
      }
      return id; 
    }

    /*!
      Return the tree's current max depth.
     */
    size_t
    max_depth() const
    {
      return max_depth_;
    }

    /*!
      Get an entity by entity id.
     */
    entity_t*
    get(
      entity_id_t id
    )
    {
      assert(id < entities_.size());
      return &(entities_[id]);
    }

    entity_t* 
    get(
        size_t id)
    {
      assert(id < entities_.size()); 
      return &(entities_[id]);
    }

    branch_t*
    get(
      branch_id_t id
    )
    { 
      auto itr = branch_map_.find(id);
      assert(itr != branch_map_.end());
      return &itr->second;
    }

    /*!
      Get the root branch (depth 0).
     */
    branch_t*
    root()
    {
      return &root_->second;
    }

    int64_t
    nonlocal_branches()
    {
      return nonlocal_branches_; 
    }

    void 
    nonlocal_branches_add()
    {
      ++nonlocal_branches_; 
    }

    /**
    * @brief Generic information for the tree topology 
    */
   friend std::ostream& operator<<(std::ostream& os,tree_topology& t )
   {
     os<<"Tree topology: "<<"#branches: "<<t.branch_map_.size()<<
       " #entities: "<<t.entities_.size();
     os <<" #root_subentities: "<<t.root()->sub_entities();
     os <<" #nonlocal_branches: "<<t.nonlocal_branches();
     return os;
   } 

  private:
    using branch_map_t = std::unordered_map<branch_id_t, branch_t,
      branch_id_hasher__<branch_int_t, dimension>>;

      branch_id_t
      to_branch_id(
        const point_t& p,
        size_t max_depth
      )
      {
        return branch_id_t(range_, p, max_depth);
      }

      branch_id_t
      to_branch_id(
        const point_t& p
      )
      {
        return branch_id_t(range_, p);
      }

      void
      insert(
        const entity_id_t& id,
        size_t max_depth
      )
      {
        auto ent = &(entities_[id]); 
        branch_id_t bid = to_branch_id(ent->coordinates(), max_depth);
        branch_t& b = find_parent(bid, max_depth);
        ent->set_branch_id_(b.id());

        b.insert(id);

        switch(b.requested_action_())
        {
          case action::none:
            break;
          case action::refine:
            refine_(b);
            break;
          default:
            assert(false && "invalid action");
        }
      }


    branch_t&
      find_parent_(
      branch_id_t bid
    )
    {
      for(;;)
      {
        auto itr = branch_map_.find(bid);
        if(itr != branch_map_.end())
        {
          return itr->second;
        }
        bid.pop();
      }
    }

    branch_t*
    find_parent(
      branch_t* b
    )
    {
      return find_parent(b.id());
    }

    branch_t*
    find_parent(
      branch_id_t bid
    )
    {
      return find_parent(bid, max_depth_);
    }

    branch_t&
    find_parent(
      branch_id_t bid,
      size_t max_depth
    )
    {
      branch_id_t pid = bid;
      pid.truncate(max_depth);

      return find_parent_(pid);
    }

     void
    refine_(
      branch_t& b
    )
    {
      // Not leaf anymore 
      branch_id_t pid = b.id();
      size_t depth = pid.depth() + 1;

      // If there are no children
      for(size_t i = 0; i < branch_t::num_children; ++i)
      {
        branch_id_t cid = pid;
        cid.push(i);
        branch_map_.emplace(cid,cid);
      }

      max_depth_ = std::max(max_depth_, depth);

      for(auto ent : b)
      {
        insert(ent, depth);
      }

      b.set_leaf(false); 
      b.clear();
      b.reset();
    }

    // helper method in coarsening
    // insert into p, coarsen all recursive children of b

    void
    coarsen_(
      branch_t* p,
      branch_t* b
    )
    {
      if(b->is_leaf())
      {
        return;
      }

      for(size_t i = 0; i < branch_t::num_children; ++i)
      {
        branch_t* ci = b->template child_<branch_t>(i);

        for(auto ent : *ci)
        {
          p->insert(ent);
          ent->set_branch_id_(p->id());
        }

        coarsen_(p, ci);
        branch_map_.erase(ci->id());
      }
    }

    void
    coarsen_(
      branch_t* p
    )
    {
      coarsen_(p, p);
      p->template into_leaf_<branch_t>();
      p->reset();
    }

    size_t
    get_queue_depth(
      thread_pool& pool
    )
    {
      size_t n = pool.num_threads();
      constexpr size_t rb = branch_int_t(1) << P::dimension;
      double bn = std::log2(double(rb));
      return std::log2(double(n))/bn + 1;
    }
 
  // Declared before, here for readability 
  //using branch_map_t = std::unordered_map<branch_id_t, branch_t,
  //    branch_id_hasher__<branch_int_t, dimension>>;

  branch_map_t branch_map_;
  size_t max_depth_;
  typename std::unordered_map<branch_id_t,branch_t,
    branch_id_hasher__<branch_int_t, dimension>>::iterator root_;
  entity_space_t entities_;
  std::array<point__<element_t, dimension>, 2> range_;
  point__<element_t, dimension> scale_;
  element_t max_scale_;
  std::vector<entity_t> entities_vector_;
  //std::vector<entity_id_t> entities_ghosts_;
  int64_t entities_vector_current_; 

  std::map<entity_id_t,entity_id_t> ghosts_id_;

  int64_t nonlocal_branches_;  
};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
