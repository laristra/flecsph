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
#include <omp.h>
#include <mpi.h>
#include <thread>

#include "flecsi/geometry/point.h"

#include "key_id.h"

#include "tree_branch.h"
#include "tree_entity.h"
#include "tree_geometry.h"
#include "entity.h"

namespace flecsi {
namespace topology {

template <class T>
bool is_unique(std::vector<T> X) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::sort(X.begin(), X.end());
  auto last = std::unique(X.begin(), X.end());
  if(last != X.end()){
    // list all
    for(auto e: X){
      std::cerr<<rank<<" "<<e<<std::endl;
    }
    //clog(trace)<<rank<<": "<<*last<<std::endl;
  }
  return last == X.end();
}

// Hasher for the branch id used in the unordered_map data structure
template<
  typename T,
  size_t D,
  class IDTYPE=key_id__<T,D>
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
  enum mpi_comm: int {
    LOCAL_REQUEST = 1,
    SOURCE_REQUEST = 2,
    MPI_DONE =3,
    SOURCE_REPLY = 5,
    MPI_RANK_DONE = 6,
    FAILED_PROBE = 7
  };


public:
  using Policy = P;

  static const size_t dimension = Policy::dimension;
  using element_t = typename Policy::element_t;
  using point_t = point__<element_t, dimension>;
  using range_t = std::array<point_t, 2>;
  using key_int_t = typename Policy::key_int_t;
  using key_t = key_id__<key_int_t,dimension>;
  using branch_id_t = key_t;
  using branch_id_vector_t = std::vector<branch_id_t>;
  using branch_t = typename Policy::branch_t;
  using branch_vector_t = std::vector<branch_t*>;
  using tree_entity_t = typename Policy::tree_entity_t;
  using entity_t = typename Policy::entity_t;
  using entity_vector_t = std::vector<tree_entity_t*>;
  using apply_function = std::function<void(branch_t&)>;
  using entity_id_vector_t = std::vector<entity_id_t>;
  using geometry_t = tree_geometry<element_t, dimension>;
  using entity_space_ptr_t = std::vector<tree_entity_t*>;

  struct filter_valid{
    bool operator()(tree_entity_t* ent) const{
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

    ghosts_entities_.resize(max_traversal);
    current_ghosts = 0;
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

    ghosts_entities_.resize(max_traversal);
    current_ghosts = 0;
  }

  /**
   * @brief Destroy the tree: empty the hash-table and destroy the entities
   * lists
   */
  ~tree_topology()
  {
    //reset_ghosts();
    branch_map_.clear();
    tree_entities_.clear();
    ghosts_id_.clear();
    ghosts_entities_.clear();
    current_ghosts = 0;
    shared_entities_.clear();
    nonlocal_branches_ = 0;
  }

  /**
  * Clean the tree topology but not the local bodies
  */
  void
  clean()
  {
    branch_map_.clear();
    tree_entities_.clear();
    ghosts_id_.clear();
    for(int i = 0 ; i <= current_ghosts; ++i)
      ghosts_entities_[i].clear();
    current_ghosts = 0;
    shared_entities_.clear();
    nonlocal_branches_ = 0;

    branch_map_.emplace(branch_id_t::root(),branch_id_t::root());
    root_ = branch_map_.find(branch_id_t::root());
    assert(root_ != branch_map_.end());
    max_depth_ = 0;
  }

  /**
  * Reset the ghosts local information for the next tree traversal
  */
  void
  reset_ghosts()
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(size == 1) return;
    //clog(trace)<<"Reset the ghosts: "<<ghosts_entities_.size()<<std::endl;
    // Remove the ghosts, all the parent have to be non local
    for(int i = 0 ; i <= current_ghosts; ++i)
    {
      for(auto& g: ghosts_entities_[i])
      {
        // Find the parent and change the status
        key_t k = g.key();
        auto& b = find_parent_(k);
        assert(!b.is_local());
        for(auto eid: b)
        {
          auto ent = get(eid);
          assert(ent->owner() != rank);
          ent->setBody(nullptr);
        }
        b.clear();
        b.set_ghosts_local(false);
        b.set_requested(false);
      }
    }
    // Same for the shared_edge
    for(auto& s: shared_entities_)
    {
      // Find the parent and change the status
      key_t k = s.key();
      auto& b = find_parent_(k);
      std::vector<entity_id_t> remove;
      for(auto eid: b)
      {
        auto ent = get(eid);
        if(ent->owner() != rank){
          ent->setBody(nullptr);
          remove.push_back(eid);
        }
      }
      for(auto r: remove){
        b.remove(r);
      }
      b.set_ghosts_local(false);
      b.set_requested(false);
      if(b.is_shared()){
        b.set_locality(branch_t::LOCAL);
      }
    }

    if(tree_entities_.size() != 0)
    {
      // Find the first ghosts
      auto start_ghosts = tree_entities_.begin();
      while(start_ghosts->is_local())
      {
        start_ghosts++;
      }
      // Remove all the associate tree branches
      if(start_ghosts != tree_entities_.end())
      tree_entities_.erase(start_ghosts,tree_entities_.end());
      // Clear ghosts informations
      ghosts_id_.clear();
      for(int i = 0 ; i <= current_ghosts; ++i)
        ghosts_entities_[i].clear();
      current_ghosts = 0;
      shared_entities_.clear();

      //mpi_tree_traversal_graphviz(11);
      share_edge();
      cofm(root(),0,false);
      //mpi_tree_traversal_graphviz(12);
    }
    //mpi_tree_traversal_graphviz(8);
    // Do the local branch sharing
  }

  /**
  * \brief Share the edge particles to my direct neighbors
  * regarding the key ordering: 0 <-> 1 <-> 2 <-> 3 for 4 processes
  */
  void
  share_edge()
  {
    // Communications 2 by 2
    // Send my highest key and my lowest key
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    // First rank%2 == 0 send to rank%2 == 1
    std::array<key_t,2> my_keys;
    std::array<key_t,2> neighbor_keys;
    size_t byte_size = sizeof(std::array<key_t,2>);

    std::vector<entity_t> received_ghosts;

    int partner = rank;
    for(int i = 0 ; i < 2; ++i)
    {
      if(rank%2 == i)
      {
        partner = rank+1;
        if(partner < size)
        {
          // Get my last key
          my_keys[0] = entities_.back().key();
          // Get the parent of this entity
          my_keys[1] = find_parent_(my_keys[0]).key();
          MPI_Send(&(my_keys[0]),byte_size,MPI_BYTE,partner,1,
            MPI_COMM_WORLD);
          MPI_Recv(&(neighbor_keys[0]),byte_size,MPI_BYTE,partner,1,
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      }else{
        partner = rank-1;
        if(partner >= 0)
        {
          // Get my last key
          my_keys[0] = entities_.front().key();
          // Get the parent of this entity
          my_keys[1] = find_parent_(my_keys[0]).key();
          MPI_Recv(&(neighbor_keys[0]),byte_size,MPI_BYTE,partner,1,
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          MPI_Send(&(my_keys[0]),byte_size,MPI_BYTE,partner,1,
            MPI_COMM_WORLD);
        }
      }
      // Entities to send to other rank
      std::vector<entity_t> edge_entities;
      if(partner >= 0 && partner < size){
        //clog(trace)<<rank<<" <-> "<<partner<<std::endl;
        // Compute the conflict
        int neighbor_parent_depth = neighbor_keys[1].depth();
        int my_parent_depth = my_keys[1].depth();
        // Handle the cases
        // The keys of the parents are the same
        if(neighbor_keys[1] == my_keys[1] )
        {
            // Send the bodies with this header
            // If this rank is the small one, go from back of vector
            int max_entities = 1<<dimension;
            if(rank < partner)
            {
              auto cur = entities_.rbegin();
              auto end = entities_.rend();
              for(;cur!=end && max_entities>0; ++cur,--max_entities)
              {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(neighbor_parent_depth);
                if(key == neighbor_keys[1])
                {
                  edge_entities.push_back(*b);
                }
              }
            // If this rank is the highest, go from the front of the vector
            }else{
              auto cur = entities_.begin();
              auto end = entities_.end();
              for(;cur!=end && max_entities>0; ++cur,--max_entities)
              {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(neighbor_parent_depth);
                if(key == neighbor_keys[1])
                {
                  edge_entities.push_back(*b);
                }
              }
            }
        }else{
          // Send nothing if the parent are not the same, other branch
        }
        // If I am the highest branch
        // In this case I just look if my body go in the neighbor branch
        if(neighbor_parent_depth < my_parent_depth)
        {
          // Is there a conflict: in this case compare the bodies
          key_t nb_key = neighbor_keys[0]; nb_key.truncate(my_parent_depth);
          key_t my_key = my_keys[0]; my_key.truncate(my_parent_depth);
          // I send the bodies that go in conflict with its body
          if(nb_key == my_key)
          {
            int max_entities = 1<<dimension;
            if(rank < partner)
            {
              auto cur = entities_.rbegin();
              auto end = entities_.rend();
              for(;cur!=end && max_entities>0; ++cur,--max_entities)
              {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(my_parent_depth);
                if(key == my_keys[1])
                {
                  edge_entities.push_back(*b);
                }
              }
            // If this rank is the highest, go from the front of the vector
            }else{
              auto cur = entities_.begin();
              auto end = entities_.end();
              for(;cur!=end && max_entities>0; ++cur,--max_entities)
              {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(my_parent_depth);
                if(key == my_keys[1])
                {
                  edge_entities.push_back(*b);
                }
              }
            }
          }else{
          }
        }
        // If my neighbor have the highest branch
        if(neighbor_parent_depth > my_parent_depth)
        {
          key_t my_key = my_keys[0]; my_key.truncate(neighbor_parent_depth);
          if(my_key == neighbor_keys[1])
          {
            // Send my branch entities
            int max_entities = 1<<dimension;
            if(rank < partner)
            {
              auto cur = entities_.rbegin();
              auto end = entities_.rend();
              for(;cur!=end && max_entities>0; ++cur,--max_entities)
              {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(neighbor_parent_depth);
                if(key == neighbor_keys[1])
                {
                  edge_entities.push_back(*b);
                }
              }
            // If this rank is the highest, go from the front of the vector
            }else{
              auto cur = entities_.begin();
              auto end = entities_.end();
              for(;cur!=end && max_entities>0; ++cur,--max_entities)
              {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(neighbor_parent_depth);
                if(key == neighbor_keys[1])
                {
                  edge_entities.push_back(*b);
                }
              }
            }
          }else{
            // Nothing to share but be prepare to receive
            // the other bodies at the right depth
            int nrefine = neighbor_parent_depth - my_parent_depth;
            int depth = my_parent_depth;
            for(int i = 0 ; i < nrefine ; ++i)
            {
              // Find my branch
              key_t branch_key = my_keys[0];
              branch_key.truncate(my_parent_depth);
              auto& b = branch_map_.find(branch_key)->second;
              refine_(b);
              ++my_parent_depth;
            }
          }
        }
      }
      // 3. Send them
      if(rank%2 == i)
      {
        partner = rank+1;
        if(partner < size)
        {
          MPI_Send(&(edge_entities[0]),sizeof(entity_t)*edge_entities.size(),
            MPI_BYTE,partner,1,MPI_COMM_WORLD);
          MPI_Status status;
          MPI_Probe(partner,1,MPI_COMM_WORLD,&status);
          // Get the count
          int nrecv = 0;
          MPI_Get_count(&status,MPI_BYTE,&nrecv);
          int offset = received_ghosts.size();
          received_ghosts.resize(offset+nrecv/sizeof(entity_t));
          MPI_Recv(&(received_ghosts[offset]),nrecv,MPI_BYTE,partner,1,
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      }else{
        partner = rank-1;
        if(partner >= 0)
        {
          MPI_Status status;
          MPI_Probe(partner,1,MPI_COMM_WORLD,&status);
          // Get the count
          int nrecv = 0;
          MPI_Get_count(&status,MPI_BYTE,&nrecv);
          int offset = received_ghosts.size();
          received_ghosts.resize(offset+nrecv/sizeof(entity_t));
          MPI_Recv(&(received_ghosts[offset]),nrecv,MPI_BYTE,partner,1,
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          MPI_Send(&(edge_entities[0]),sizeof(entity_t)*edge_entities.size(),
            MPI_BYTE,partner,1,MPI_COMM_WORLD);
        }
      }
      //clog(trace)<<rank<<" received "<<received_ghosts.size()<<std::endl;
    }

    for(auto g: received_ghosts){
      shared_entities_.push_back(g);
    }
    for(auto& g : shared_entities_)
    {
      //clog(trace)<<"Inserting in the tree "<<g.key()<<std::endl;
      auto* bi = &(g);
      assert(bi->mass()!=0.);
      assert(bi->key().value_() != 0);
      if(bi->key() > entities().back().key()){
        partner = rank +1;
      }else{
        partner = rank -1;
      }
      auto id = make_entity(bi->key(),bi->coordinates(),nullptr,partner,
        bi->mass(),bi->id(),bi->radius());
      insert(id);
      get(id)->setBody(bi);
      // Set the ghosts local in this case
    }
    for(auto& g : shared_entities_)
    {
      find_parent_(g.key()).set_ghosts_local(true);
    }
  }

  /**
  * \brief Change the range of the tree topology
  */
  void
  set_range(
    const range_t& range)
  {
    range_ = range;
    max_scale_ = range_[1][0] - range_[0][0];

    for(size_t d = 0; d < dimension; ++d)
    {
      scale_[d] = range_[1][d] - range_[0][d];
      max_scale_ = std::max(max_scale_, scale_[d]);
    }
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
    Return an index space containing all non-removed entities.
   */
  std::vector<tree_entity_t>&
  tree_entities()
  {
    return tree_entities_;
  }

  std::vector<entity_t>&
  entities()
  {
    return entities_;
  }

  template<
    typename E>
  entity_t&
  entity(E e){
    return entities_[static_cast<int>(e)];
  }

  //std::vector<entity_id_t>&
  //tree_entities_ghosts()
  //{
  //  return tree_entities_ghosts_;
  //}

  tree_entity_t*
  get_ghost(
      const size_t& global_id)
  {
    auto local_id_itr = ghosts_id_.find(global_id);
    assert(local_id_itr != ghosts_id_.end());
    auto local_id = local_id_itr->second;
    return &(tree_entities_[local_id]);
  }

  /*!
    Insert an entity into the lowest possible branch division.
   */
  void
  insert(
    const entity_id_t& id
  )
  {
    //std::cout<<"Inserting -------------> "<<id<<std::endl;
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
          if(next == nullptr) continue;
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
          if(next == nullptr) continue;
          search_list.push_back(next);
          stk.push(next);
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
  get_sub_entities(
    branch_t * start,
    std::vector<entity_t>& search_list)
  {
    std::stack<branch_t*> stk;
    stk.push(start);

    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        for(auto id: *c){
          auto child = this->get(id);
          if(child->getBody() != nullptr){
            search_list.push_back(*(child->getBody()));
          }
        }
      }else{
        for(int i=0; i<(1<<dimension);++i){
          branch_t * next = child(c,i);
          if(next == nullptr) continue;
          stk.push(next);
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
      if(c->is_leaf()){
        search_list.push_back(c);
      }else{
        if(c->sub_entities() < criterion){
          search_list.push_back(c);
        }else{
          for(int i=0; i<(1<<dimension);++i){
            branch_t * next = child(c,i);
            if(next == nullptr) continue;
            stk.push(next);
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
            if(next == nullptr) continue;
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
  apply_sub_cells(
      branch_t * b,
      int64_t ncritical,
      bool do_square,
      EF&& ef,
      ARGS&&... args)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    std::stack<branch_t*> stk;
    stk.push(b);

    std::vector<branch_t*> work_branch;

    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf() && c->is_local()){
        work_branch.push_back(c);
      }else{
        if(c->is_local() && (int64_t)c->sub_entities() < ncritical){
          work_branch.push_back(c);
        }else{
          for(int i=0; i<(1<<dimension);++i){
            branch_t * next = child(c,i);
            if(next == nullptr) continue;
            if(next->is_local()){
              stk.push(next);
            }
          }
        }
      }
    }

    std::vector<branch_t*> remaining_branches;
    int done = omp_get_max_threads();
    #pragma omp parallel num_threads(omp_get_max_threads()+1)
    {
      int nthreads = omp_get_num_threads();
      int tid = omp_get_thread_num();

      if(tid == nthreads-1){
        handle_requests();
      }else{
        traverse_branches(do_square,work_branch,remaining_branches,false,
          ef,std::forward<ARGS>(args)...);
        // Wait for other threads
        #pragma omp critical
        {
          if(--done == 0)
            MPI_Send(NULL, 0, MPI_INT, rank, MPI_DONE,MPI_COMM_WORLD);
        }
      }
    }
    //clog(trace)<<"Merged done"<<std::endl<<std::flush;

    // Check if no message remainig
#ifdef DEBUG
    int flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
      MPI_STATUS_IGNORE);
    assert(flag == 0);
#endif
    for(size_t i = 0; i < ghosts_entities_[current_ghosts].size(); ++i)
    {
      entity_t& g = ghosts_entities_[current_ghosts][i];
      auto id = make_entity(g.key(),g.coordinates(),
        nullptr,g.owner(),g.mass(),g.id(),g.radius());
      // Assert the parent exists and is non local
      assert(!find_parent_(g.key()).is_local());
      insert(id);
      auto nbi = get(id);
      nbi->setBody(&g);
      assert(nbi->global_id() == g.id());
      assert(nbi->getBody() != nullptr);
      // Set the parent to local for the search
      find_parent_(g.key()).set_ghosts_local(true);
    }
    ++current_ghosts;
    assert(current_ghosts < max_traversal);
    cofm(root(), 0, false);

    std::vector<branch_t*> ignore;
    //clog(trace)<<"Cmopute the remaining branches"<<std::endl<<std::flush;
    #pragma omp parallel
    {
      traverse_branches(do_square,remaining_branches,ignore,true,
        ef,std::forward<ARGS>(args)...);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  } // apply_sub_cells

  template<
    typename EF,
    typename... ARGS
  >
  void
  traverse_branches(
    bool do_square,
    std::vector<branch_t*>& work_branch,
    std::vector<branch_t*>& non_local_branches,
    bool assert_local,
    EF&& ef,
    ARGS&&... args)
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int nelem = work_branch.size();
    int td = omp_get_thread_num();
    int nt = omp_get_num_threads();
    if(!assert_local){
      --nt;
    }
    int start = td*nelem/nt;
    int end = nelem*(td+1)/nt;
    if(td == nt-1){
      end = nelem;
    }
    for(size_t i = start; i < end; ++i){
      std::vector<branch_t*> inter_list;
      std::vector<branch_t*> requests_branches;
      if(sub_cells_inter(work_branch[i],inter_list,requests_branches)){
        force_calc(work_branch[i],inter_list,do_square,
               ef,std::forward<ARGS>(args)...);
      }else{
        assert(!assert_local);
        std::vector<key_t> send;
        #pragma omp critical
        {
          non_local_branches.push_back(work_branch[i]);
          // Send branch key to request handler
          for(auto b: requests_branches){
            assert(b->owner() < size && b->owner() >= 0);
            assert(b->owner() != rank);
            if(!b->requested()){
              send.push_back(b->id());
              b->set_requested(true);
            }
          }
        }
        MPI_Request request;
        MPI_Send(&(send[0]),send.size()*sizeof(key_t),MPI_BYTE,rank,
          LOCAL_REQUEST,MPI_COMM_WORLD);
      }
    }
  }


  template<
    typename FC,
    typename DFCDR,
    typename DFCDRDR,
    typename C2P
  >
  void
  apply_fmm(
      branch_t * b,
      double maxmasscell,
      const double MAC,
      FC f_fc,
      DFCDR f_dfcdr,
      DFCDRDR f_dfcdrdr,
      C2P f_c2p
    )
  {

    MPI_Barrier(MPI_COMM_WORLD);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    rank || clog(trace) << "FMM : maxmasscell: "<<maxmasscell<<" MAC: "<<
      MAC<<std::endl;

    std::stack<branch_t*> stk;
    stk.push(b);
    std::vector<branch_t*> work_branch;
    // Only work on local branches
    // 1. Gather the cells
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf() && c->is_local()){
        work_branch.push_back(c);
      }else{
        if(c->is_local() && c->mass() <= maxmasscell){
          work_branch.push_back(c);
        }else{
          for(int i=0; i<(1<<dimension);++i){
            branch_t * next = child(c,i);
            if(next == nullptr) continue;
            if(next->is_local()){
              stk.push(next);
            }
          }
        }
      }
    }

    std::vector<branch_t*> remaining_branches;
    int done = omp_get_max_threads();
    #pragma omp parallel num_threads(omp_get_max_threads()+1)
    {
      int nthreads = omp_get_num_threads();
      int tid = omp_get_thread_num();

      if(tid == nthreads-1){
        handle_requests();
      }else{
        traversal_fmm(work_branch,remaining_branches,MAC,
          f_fc,f_dfcdr,f_dfcdrdr,f_c2p);
        // Wait for other threads
        #pragma omp critical
        {
          if(--done == 0){
            MPI_Send(NULL, 0, MPI_INT, rank, MPI_DONE,MPI_COMM_WORLD);
          }
          assert(done >= 0);
        }
      }
    }

    // Check if no message remainig
#ifdef DEBUG
    int flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
      MPI_STATUS_IGNORE);
    assert(flag == 0);
#endif
    for(size_t i = 0; i < ghosts_entities_[current_ghosts].size(); ++i)
    {
      entity_t& g = ghosts_entities_[current_ghosts][i];
      auto id = make_entity(g.key(),g.coordinates(),
        nullptr,g.owner(),g.mass(),g.id(),g.radius());
      // Assert the parent exists and is non local
      assert(!find_parent_(g.key()).is_local());
      insert(id);
      auto nbi = get(id);
      nbi->setBody(&g);
      assert(nbi->global_id() == g.id());
      assert(nbi->getBody() != nullptr);
      // Set the parent to local for the search
      find_parent_(g.key()).set_ghosts_local(true);
    }
    ++current_ghosts;
    assert(current_ghosts < max_traversal);
    cofm(root(), 0, false);

    // compute the particle to particle interaction on each sub-branches
    #pragma omp parallel for
    for(int i = 0 ; i < work_branch.size(); ++i)
    {
      traversal_p2p(work_branch[i],MAC,f_fc);
    }
  } // apply_sub_cells

  template<
    typename FC,
    typename DFCDR,
    typename DFCDRDR,
    typename C2P
  >
  void
  traversal_fmm(
    std::vector<branch_t*>& work_branch,
    std::vector<branch_t*>& non_local_branches,
    const double MAC,
    FC f_fc,
    DFCDR f_dfcdr,
    DFCDRDR f_dfcdrdr,
    C2P f_c2p
  )
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int nelem = work_branch.size();
    int td = omp_get_thread_num();
    int nt = omp_get_num_threads();
    // Do not count the thread handling the requests
    --nt;
    int start = td*nelem/nt;
    int end = nelem*(td+1)/nt;
    if(td == nt-1){
      end = nelem;
    }
    for(size_t i = start; i < end; ++i){
      std::vector<branch_t*> inter_list;
      std::vector<branch_t*> requests_branches;
      // Sub traversal and compute the MAC

      if(traversal_c2c_c2p(work_branch[i],requests_branches,MAC,
        f_fc,f_dfcdr,f_dfcdrdr,f_c2p))
      {
        std::vector<key_t> send;
        #pragma omp critical
        {
          non_local_branches.push_back(work_branch[i]);
          // Send branch key to request handler
          for(auto b: requests_branches){
            assert(b->owner() < size && b->owner() >= 0);
            assert(b->owner() != rank);
            if(!b->requested()){
              send.push_back(b->id());
              b->set_requested(true);
            }
          }
        }
        MPI_Request request;
        MPI_Send(&(send[0]),send.size()*sizeof(key_t),MPI_BYTE,rank,
          LOCAL_REQUEST,MPI_COMM_WORLD);
      }
    }
  }

  template<
    typename FC,
    typename DFCDR,
    typename DFCDRDR,
    typename F_C2P
  >
  bool
  traversal_c2c_c2p(
    branch_t* b,
    std::vector<branch_t*>& non_local,
    const double MAC,
    FC f_fc, DFCDR f_dfcdr,DFCDRDR f_dfcdrdr, F_C2P f_c2p
  )
  {
    point_t fc = 0.;
    double dfcdr[9] = {0.};
    double dfcdrdr[27] = {0.};

    std::stack<branch_t*> stk;
    stk.push(root());
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        if(geometry_t::box_MAC(
          b->coordinates(),c->coordinates(),
          c->bmin(),c->bmax(),MAC)
        ) {
          // Compute the matrices
          f_fc(fc,b->coordinates(),c->coordinates(),c->mass());
          f_dfcdr(dfcdr,b->coordinates(),c->coordinates(),c->mass());
          f_dfcdrdr(dfcdrdr,b->coordinates(),c->coordinates(),c->mass());
        }else if(!c->is_local() && !c->ghosts_local()){
          // Request the sub-particles of this leaf
          non_local.push_back(c);
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(c,i);
          if(branch == nullptr) continue;
          if(geometry_t::box_MAC(
            b->coordinates(),
            branch->coordinates(),
            branch->bmin(),
            branch->bmax(),MAC))
          {
            stk.push(branch);
          }
        }
      }
    }

    // Propagate this information to the sub-particles for C2P
    stk.push(b);
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        // For every bodies propagate
        for(auto id: *b){
          auto child = &(tree_entities_[id]);
          if(child->is_local())
            f_c2p(fc,dfcdr,dfcdrdr,b->coordinates(),child);
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(c,i);
          if(branch == nullptr) continue;
          stk.push(branch);
        }
      }
    }
    return non_local.size() > 0;
  }

  template<typename F_FC>
  void
  traversal_p2p(
    branch_t* b,
    element_t MAC,
    F_FC f_fc
  )
  {
    // Collect all sub-entities
    std::vector<tree_entity_t*> local_entities;
    std::stack<branch_t*> stk;
    stk.push(b);
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf() && c->is_local()){
        for(auto id: *b){
          auto child = &(tree_entities_[id]);
          if(child->is_local()){
            local_entities.push_back(child);
          }
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(c,i);
          if(branch == nullptr) continue;
          stk.push(branch);
        }
      }
    }

    stk.push(root());
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        if(!geometry_t::box_MAC(
          b->coordinates(),c->coordinates(),
          c->bmin(),c->bmax(),MAC)
        )
        {
          for(auto id: *b){
            auto child = &(tree_entities_[id]);
            for(auto& loc: local_entities)
            {
              // N square computation
              point_t fc;
              loc->getBody()->setAcceleration(
                loc->getBody()->getAcceleration() +
                f_fc(fc,loc->coordinates(),child->coordinates(),child->mass())
              );
            }
          }
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(c,i);
          if(branch == nullptr) continue;
          if(!geometry_t::box_MAC(
            b->coordinates(),
            branch->coordinates(),
            branch->bmin(),
            branch->bmax(),MAC))
          {
            stk.push(branch);
          }
        }
      }
    }
  }

/**
*
*/
  void
  handle_requests()
  {
    int max_size = 500;
    int max_requests = 100;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    bool done_traversal = false;
    bool done = false;

    // Maintain a request array for all neighbors
    std::vector<std::vector<std::vector<key_t>>> requests(size);
    std::vector<std::vector<std::vector<entity_t>>> reply(size);
    std::vector<int> current_requests(size,0);
    std::vector<int> current_reply(size,0);
    for(int i = 0 ; i < size; ++i){
      requests[i].resize(max_requests);
      reply[i].resize(max_requests);
    }

    std::vector<MPI_Request> mpi_requests(max_requests);
    int current_mpi_requests = 0;
    std::vector<MPI_Request> mpi_replies(max_requests);
    int current_mpi_replies = 0;
    int request_counter = 0;
    std::vector<bool> rank_done(size,false);
    bool done_rank = false;

    while(!done)
    {
      MPI_Status status;
      // Wait on probe

      int source, tag, nrecv;
      if(!done_traversal){
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        source = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        nrecv = 0;
        MPI_Get_count(&status, MPI_BYTE, &nrecv);
      }else{
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
            &status);
        if(flag){
          source = status.MPI_SOURCE;
          tag = status.MPI_TAG;
          nrecv = 0;
          MPI_Get_count(&status, MPI_BYTE, &nrecv);
        }else{
          tag = FAILED_PROBE;
        }
      }
      //if(tag != LOCAL_REQUEST)
      //  clog(trace)<<rank<<" received: "<<tag<<" from "<< source <<std::endl;
      switch(tag)
      {
        case MPI_RANK_DONE:
        {
          MPI_Recv(NULL,0,MPI_INT,source,MPI_RANK_DONE,MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
          rank_done[source] = true;
        }
        break;
        // ------------------------------------------------------------------ //
        //     Another rank replied to my entities request                    //
        // ------------------------------------------------------------------ //
        case SOURCE_REPLY:
        {
          assert(nrecv != 0);
          std::vector<entity_t> received(nrecv/sizeof(entity_t));
          MPI_Recv(&(received[0]),nrecv,MPI_BYTE,source,
            SOURCE_REPLY,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          //assert(is_unique(received));
          ghosts_entities_[current_ghosts].insert(
            ghosts_entities_[current_ghosts].end(),received.begin(),
              received.end());
          //assert(is_unique(ghosts_entities_[current_ghosts]));
          --request_counter;
        }
        break;
        // ------------------------------------------------------------------ //
        //    Another rank request information for entities                   //
        // ------------------------------------------------------------------ //
        case SOURCE_REQUEST:
        {
          std::vector<key_t> received(nrecv/sizeof(key_t));
          MPI_Recv(&(received[0]), nrecv, MPI_BYTE, source, SOURCE_REQUEST,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          //assert(is_unique(received));
          for(auto k: received){
            auto branch = branch_map_.find(k);
            assert(branch != branch_map_.end());
            get_sub_entities(&(branch->second),
              reply[source][current_reply[source]]);
          }
          //assert(is_unique(reply[source][current_reply[source]]));
          int ncount = reply[source][current_reply[source]].size();
          MPI_Isend(&(reply[source][current_reply[source]++][0]),
            ncount*sizeof(entity_t),MPI_BYTE,source,SOURCE_REPLY,MPI_COMM_WORLD,
            &(mpi_replies[current_mpi_replies++]));
          if(current_mpi_replies > max_requests){
              clog(error)<<rank<<
                ": Exceeding number of replies requests"<<std::endl;
          }
        }
        break;
        // ------------------------------------------------------------------ //
        //        A local thread requested a distant particles                //
        // ------------------------------------------------------------------ //
        case LOCAL_REQUEST:
        {
          std::vector<key_t> keys(nrecv/sizeof(key_t));
          MPI_Recv(&(keys[0]), nrecv, MPI_BYTE, rank, LOCAL_REQUEST,
            MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          for(auto k: keys){
            auto itr = branch_map_.find(k); assert(itr != branch_map_.end());
            branch_t* branch = &(itr->second); assert(branch->requested());
            int owner = branch->owner(); assert(owner != rank);
            // Add this request to vector
            requests[owner][current_requests[owner]].push_back(k);
            if(requests[owner][current_requests[owner]].size() >= max_size){
              //assert(is_unique(requests[owner][current_requests[owner]]));
              MPI_Isend(&(requests[owner][current_requests[owner]++][0]),
                max_size*sizeof(key_t),MPI_BYTE,owner,SOURCE_REQUEST,
                MPI_COMM_WORLD,&(mpi_requests[current_mpi_requests++]));
              if(current_mpi_requests > max_requests){
                  clog(error)<<rank<<
                    ": Exceeding number of requests requests"<<std::endl;
              }
              ++request_counter;
            }
          }
        }
        break;
        // ------------------------------------------------------------------ //
        //        The OpenMP threads are done, send the last requests         //
        // ------------------------------------------------------------------ //
        case MPI_DONE:
        {
          // First time, send remaining requests
          MPI_Recv(NULL,0,MPI_INT,source,MPI_DONE,MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
          for(int i = 0 ; i < size ; ++i){
            if(i == rank ) continue;
            if(requests[i][current_requests[i]].size() > 0)
            {
              int nsend = requests[i][current_requests[i]].size();
              MPI_Isend(&(requests[i][current_requests[i]][0]),
                nsend*sizeof(key_t),MPI_BYTE,i,SOURCE_REQUEST,MPI_COMM_WORLD,
                &(mpi_requests[current_mpi_requests++]));
              if(current_mpi_requests > max_requests)
              {
                  clog(error)<<rank<<
                    ": Exceeding number of requests requests"<<std::endl;
              }
              ++request_counter;
            }
          }
          done_traversal = true;
        }
        break;

        case FAILED_PROBE:
          break;

        default:
          {assert(false);}
        break;
      }

      if(done_traversal)
      {
        // Check if all requests have been answered
        done = request_counter==0?true:false;
        bool done_requests = true;
        int flag;
        for(size_t i = 0 ; i < current_mpi_requests; ++i){
          MPI_Test(&(mpi_requests[i]),&flag,MPI_STATUS_IGNORE);
          done_requests = done_requests && flag;
        }
        done = done && done_requests;
        bool done_replies = true;
        for(size_t i = 0 ; i < current_mpi_replies; ++i){
          MPI_Test(&(mpi_replies[i]),&flag,MPI_STATUS_IGNORE);
          done_replies = done_requests && flag;
        }
        done = done && done_replies;
        if(done && (!done_rank))
        {
          done_rank = true;
          // Send rank done
          for(int i = 0 ; i < size; ++i){
            if(i == rank) continue;
            MPI_Isend(NULL,0,MPI_INT,i,MPI_RANK_DONE,MPI_COMM_WORLD,
                &(mpi_requests[current_mpi_requests++]));
            if(current_mpi_requests > max_requests)
            {
                clog(error)<<rank<<
                  ": Exceeding number of requests requests"<<std::endl;
            }
          }
          rank_done[rank] = true;
        }
        for(size_t i = 0 ; i < size; ++i){
          done = done && rank_done[i];
        }
      }
    }
    //clog(trace)<<"Handler done"<<std::endl;
  }

  void
  cofm(branch_t * start, element_t epsilon, bool local)
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
          if(next == nullptr) continue;
          stk1.push(next);
        }
      }
    }
    // Unstack stack 2
    while(!stk2.empty())
    {
      branch_t * cur = stk2.top();
      stk2.pop();
      update_COM(cur,epsilon,local);
    }
  }

   // Functions for the tree traversal
   void
   update_COM(
        branch_t * b,
        element_t epsilon = element_t(0),
        bool local_only = false)
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);

      //typename branch_t::b_locality locality = branch_t::NONLOCAL;
      element_t mass = element_t(0);
      point_t bmax{};
      point_t bmin{};
      element_t radius = 0.;
      point_t coordinates = point_t{};
      uint64_t nchildren = 0;
      int owner = b->owner();
      for(size_t d = 0 ; d < dimension ; ++d){
        bmax[d] = -DBL_MAX;
        bmin[d] = DBL_MAX;
      }
      bool full_nonlocal = true;
      bool full_local = true;
      if(b->is_leaf()){
        // For local branches, compute the radius
        if(b->is_local()){
          for(auto child: *b)
          {
            auto ent = get(child);
            owner = ent->owner();
            if(local_only && !ent->is_local()){
              continue;
            }
            if(ent->owner() == rank){
              full_nonlocal = false;
            }else{
              full_local = false;
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
            auto ent = get(child);
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
          nchildren = b->sub_entities();
        }
        // Locality for leaves
        if(full_nonlocal && !full_local){
          b->set_owner(owner);
        }
        if(b->owner() == rank && full_local){
          b->set_locality(branch_t::LOCAL);
        }else if(b->owner() == rank && !full_local){
          b->set_locality(branch_t::SHARED);
        }else{
          b->set_locality(branch_t::NONLOCAL);
        }
      }else{
        bool local = false;
        bool nonlocal = false;
        for(int i = 0 ; i < (1<<dimension); ++i)
        {
          auto branch = child(b,i);
          if(branch == nullptr) continue;

          nchildren+=branch->sub_entities();
          mass += branch->mass();
          if(branch->locality() == branch_t::LOCAL){
            local = true;
          }
          if(branch->locality() == branch_t::NONLOCAL){
            nonlocal = true;
          }
          if(branch->locality() == branch_t::SHARED){
            local = nonlocal = true;
          }
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
          auto branch = child(b,i);
          if(branch == nullptr) continue;

          radius = std::max(
              radius,
              distance(coordinates,branch->coordinates()) + branch->radius());
        }
        if(local && nonlocal){
          b->set_locality(branch_t::SHARED);
        }
        if(local && !nonlocal){
          b->set_locality(branch_t::LOCAL);
        }
        if(!local && nonlocal){
          b->set_locality(branch_t::NONLOCAL);
        }
      }
      b->set_radius(radius);
      b->set_sub_entities(nchildren);
      b->set_coordinates(coordinates);
      b->set_mass(mass);
      b->set_bmin(bmin);
      b->set_bmax(bmax);
      if(nchildren == 0){
        b->set_locality(branch_t::EMPTY);
      }
      if(!b->is_local()) nonlocal_branches_add();
    }


  bool
  sub_cells_inter(
    branch_t* b,
    std::vector<branch_t*>& inter_list,
    std::vector<branch_t*>& non_local)
  {
    std::stack<branch_t*> stk;
    stk.push(root());
    while(!stk.empty()){
      branch_t* c = stk.top();
      stk.pop();
      if(c->is_leaf()){
        if(c->is_local() || c->ghosts_local()){
          inter_list.push_back(c);
        }else{
          non_local.push_back(c);
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(c,i);
          if(branch == nullptr) continue;
          if(geometry_t::intersects_box_box(
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
    return (non_local.size() == 0);
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
          if(next == nullptr) continue;
          stk.push(next);
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
      tree_entity_t* ent,
      std::vector<branch_t*>& inter_list,
      EF&& ef,
      ARGS&&... args)
  {
    std::vector<tree_entity_t*> neighbors;
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
      tree_entity_t* ent,
      std::vector<branch_t*>& inter_list,
      EF&& ef,
      ARGS&&... args)
  {
    std::vector<tree_entity_t*> neighbors;
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
          auto child = &(tree_entities_[id]);
          // Check if in radius
          if(ef(center,child->coordinates(),radius,child->h())){
            ents.push_back(child);
          }
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(b,i);
          if(branch == nullptr ) continue;
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
          auto child = &(tree_entities_[id]);
          // Check if in box
          if(ef(min,max,child->coordinates(),child->h())){
            ents.push_back(child);
          }
        }
      }else{
        for(int i=0 ; i<(1<<dimension);++i){
          auto branch = child(b,i);
          if(branch == nullptr) continue;
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
      tree_entities_.emplace_back(std::forward<Args>(args)...);
      auto ent = &(tree_entities_.back());
      // Size -1 to start at 0
      entity_id_t id = tree_entities_.size()-1;
      ent->set_id_(id);
      if(!ent->is_local()){
        ghosts_id_.insert(std::make_pair(ent->global_id(),id));
      }
      return id;
    }

    /**
    * Insert directly a branch (certainly remote) in the tree
    */
    void
    insert_branch(
      const point_t& coordinates,
      const element_t& mass,
      const point_t& bmin,
      const point_t& bmax,
      const key_t& key,
      const int& owner,
      const size_t& sub_entities
    ){
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      assert(rank != owner);
      //clog(trace)<<rank<<" inserting branch: "<<key<<std::endl;
      // Check if this key already exists
      auto itr = branch_map_.find(key);
      // Case 1, branch does not exists localy
      if(itr == branch_map_.end()){
        // Add the missing parents
        key_t pk = key;
        pk.pop();
        while(branch_map_.find(pk) == branch_map_.end()){
          //clog(trace)<<rank<<" adding parent: "<<pk<<std::endl;
          branch_map_.emplace(pk,pk);
          itr = branch_map_.find(pk);
          itr->second.set_ghosts_local(false);
          itr->second.set_coordinates(coordinates);
          itr->second.set_mass(mass);
          itr->second.set_bmin(bmin);
          itr->second.set_bmax(bmax);
          itr->second.set_owner(owner);
          itr->second.set_sub_entities(sub_entities);
          itr->second.set_locality(branch_t::NONLOCAL);
          itr->second.set_leaf(false);
          pk.pop();
        }
        // Set upper level not to leave
        branch_map_.find(pk)->second.set_leaf(false);

        branch_map_.emplace(key,key);
        itr = branch_map_.find(key);
        itr->second.set_ghosts_local(false);
        itr->second.set_coordinates(coordinates);
        itr->second.set_mass(mass);
        itr->second.set_bmin(bmin);
        itr->second.set_bmax(bmax);
        itr->second.set_owner(owner);
        itr->second.set_sub_entities(sub_entities);
        itr->second.set_locality(branch_t::NONLOCAL);
        itr->second.set_leaf(true);
      }else{

        //clog(trace) <<rank << " key already exists" << std::endl;
        if(itr->second.sub_entities() == 0){
          //clog(trace) <<rank<< " empty branch" << coordinates<< std::endl;
          itr->second.set_ghosts_local(false);
          itr->second.set_coordinates(coordinates);
          itr->second.set_mass(mass);
          itr->second.set_bmin(bmin);
          itr->second.set_bmax(bmax);
          itr->second.set_owner(owner);
          itr->second.set_sub_entities(sub_entities);
          itr->second.set_locality(branch_t::NONLOCAL);
          itr->second.set_leaf(true);
        }else{
          if(itr->second.owner() == rank)
            assert(itr->second.is_shared());
          else
            assert(!itr->second.is_shared());
          // Check if same children exists
        }
      }
      // Add this branch if does not exists
    }

    void
    compute_keys()
    {
      #pragma omp parallel for
      for(size_t i = 0; i < entities_.size(); ++i){
        entities_[i].set_key(
            key_t(range_,entities_[i].coordinates()));
      }
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
    tree_entity_t*
    get(
      entity_id_t id
    )
    {
      assert(id < tree_entities_.size());
      return &(tree_entities_[id]);
    }

    tree_entity_t*
    get(
        size_t id)
    {
      assert(id < tree_entities_.size());
      return &(tree_entities_[id]);
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


    std::vector<entity_t>&
    ghosts_entities()
    {
      return ghosts_entities_[current_ghosts];
    }

    std::vector<entity_t>&
    shared_entities()
    {
      return shared_entities_;
    }

    /**
    * @brief Generic information for the tree topology
    */
   friend std::ostream& operator<<(std::ostream& os,tree_topology& t )
   {
     os<<"Tree topology: "<<"#branches: "<<t.branch_map_.size()<<
       " #entities: "<<t.tree_entities_.size();
     os <<" #root_subentities: "<<t.root()->sub_entities();
     os <<" #nonlocal_branches: "<<t.nonlocal_branches();
     return os;
   }


   /**
    * @brief      Export to a file the current tree in memory
    * This is useful for small number of particles to help representing the tree
    *
    * @param      tree   The tree to output
    * @param      range  The range of the particles, use to construct entity_key
    */
   void
   mpi_tree_traversal_graphviz(
     int num)
   {
     int rank = 0;
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
     rank || clog(trace)<<rank<<" outputing tree file #"<<num<<std::endl;

     char fname[64];
     sprintf(fname,"output_graphviz_%02d_%02d.gv",rank,num);
     std::ofstream output;
     output.open(fname);
     output<<"digraph G {"<<std::endl<<"forcelabels=true;"<<std::endl;

     // Add the legend
     output<<"branch [label=\"branch\" xlabel=\"sub_entities,owner,requested,"
      "ghosts_local\"]"<<std::endl;

     std::stack<branch_t*> stk;
     // Get root
     auto rt = root();
     stk.push(rt);

     while(!stk.empty()){
       branch_t* cur = stk.top();
       stk.pop();
       if(!cur->is_leaf()){
           output<<cur->id()<<" [label=\""<<cur->id()<< "\", xlabel=\"" <<
           cur->sub_entities()<<" - "<< cur->owner() <<" - "<< cur->requested()
           <<" - "<<cur->ghosts_local()<<"\"];"<<std::endl;
         switch (cur->locality())
         {
           case 1:
             output<<cur->id()<<" [shape=circle,color=blue]"<<std::endl;
             break;
           case 2:
             output<<cur->id()<<" [shape=circle,color=red]"<<std::endl;
             break;
           case 3:
             output<<cur->id()<<" [shape=circle,color=green]"<<std::endl;
             break;
           default:
             output<<cur->id()<<" [shape=circle,color=black]"<<std::endl;
             break;
         }

         // Add the child to the stack and add for display
         for(size_t i=0;i<(1<<dimension);++i)
         {
           auto br = child(cur,i);
           if(br == nullptr) continue;
           stk.push(br);
             output<<std::oct<<cur->id()<<"->"<<br->id()<<std::dec<<std::endl;
         }
       }else{
         output<<cur->id()<<" [label=\""<<cur->id()<< "\", xlabel=\"" <<
           cur->sub_entities()<< " - "<< cur->owner() <<" - "<< cur->requested()
           <<" - "<<cur->ghosts_local()<<"\"];"<<std::endl;
         switch (cur->locality())
         {
           case 1:
             output<<cur->id()<<" [shape=circle,color=blue]"<<std::endl;
             break;
           case 2:
             output<<cur->id()<<" [shape=circle,color=red]"<<std::endl;
             break;
           case 3:
             output<<cur->id()<<" [shape=circle,color=green]"<<std::endl;
             break;
           default:
             output<<cur->id()<<" [shape=circle,color=black]"<<std::endl;
             break;
         }
         for(auto ent: *cur)
         {
           auto e = get(ent);
           key_t key(range(),e->coordinates());
           key.truncate(max_depth()+2);

           output<<key<<" [label=\""<<key << "\", xlabel=\"" <<
            e->owner() <<" - "<<e->global_id()<<"\"];"<<std::endl;

           output<<cur->id()<<"->"<<key<<std::endl;
           switch (e->locality())
           {
             case 2:
               output<<key<<" [shape=box,color=green]"<<std::endl;
               break;
             case 3:
               output<<key<<" [shape=box,color=black]"<<std::endl;
               break;
             case 1:
               output<<key<<" [shape=box,color=red]"<<std::endl;
               break;
             default:
               output<<key<<" [shape=circle,color=red]"<<std::endl;
               break;
           }
           output<<std::dec;
         }
       }
     }
     output<<"}"<<std::endl;
     output.close();
     //MPI_Barrier(MPI_COMM_WORLD);
   }

  private:
    using branch_map_t = std::unordered_map<branch_id_t, branch_t,
      branch_id_hasher__<key_int_t, dimension>>;

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

        // Find parent of the id
        auto ent = &(tree_entities_[id]);
        branch_id_t bid = ent->get_entity_key();
        branch_t& b = find_parent(bid, max_depth);
        // Check if it is a leaf
        if(!b.is_leaf())
        {
          //clog(trace)<<"Not leaf"<<std::endl;
          // Create the branch
          int depth = b.id().depth()+1;
          bid.truncate(depth);
          // Insert this branch and reinsert
          branch_map_.emplace(bid,bid);
          //clog(trace)<<"Creating sub parent: "<<bid<<std::endl;
          branch_map_.find(bid)->second.set_leaf(true);
          branch_map_.find(bid)->second.insert(id);
        }else{
          // Check if there is a conflict with existing children
          bool conflict = false;
          int e_depth = b.id().depth()+1;
          key_t new_key = ent->key();
          new_key.truncate(e_depth);
          for(auto id: b)
          {
            auto e = get(id);
            key_t k = e->key();
            k.truncate(e_depth);
            if(k == new_key)
            {
              conflict = true;
            }
          }
          // Conflict with a children
          if(conflict){
            refine_(b);
            insert(id);
          }else{
            b.insert(id);
          }
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
      branch_id_t pid = b.id();
      size_t depth = pid.depth() + 1;

      // For every children
      for(auto ent: b)
      {
        key_t k = get(ent)->key();
        k.truncate(depth);
        branch_map_.emplace(k,k);
        insert(ent,depth);
      }
      max_depth_ = std::max(max_depth_, depth);

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
      constexpr size_t rb = key_int_t(1) << P::dimension;
      double bn = std::log2(double(rb));
      return std::log2(double(n))/bn + 1;
    }

  // Declared before, here for readability
  //using branch_map_t = std::unordered_map<branch_id_t, branch_t,
  //    branch_id_hasher__<key_int_t, dimension>>;

  branch_map_t branch_map_;
  size_t max_depth_;
  typename std::unordered_map<branch_id_t,branch_t,
      branch_id_hasher__<key_int_t, dimension>>::iterator root_;
  range_t range_;
  point__<element_t, dimension> scale_;
  element_t max_scale_;
  std::map<entity_id_t,entity_id_t> ghosts_id_;
  std::vector<tree_entity_t> tree_entities_;
  std::vector<entity_t> entities_;

  const size_t max_traversal = 5;
  std::vector<std::vector<entity_t>> ghosts_entities_;
  size_t current_ghosts = 0;

  std::vector<entity_t> shared_entities_;

  int64_t nonlocal_branches_;
};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
