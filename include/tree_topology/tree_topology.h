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

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <float.h>
#include <functional>
#include <iostream>
#include <map>
#include <math.h>
#include <mpi.h>
#include <mutex>
#include <omp.h>
#include <set>
#include <stack>
#include <thread>
#include <unordered_map>
#include <vector>

#include "flecsi/data/data_client.h"
#include "flecsi/geometry/point.h"

#include "tree_geometry.h"
#include "tree_types.h"


namespace flecsi {
namespace topology {

/*!
  The tree topology is parameterized on a policy P which defines its branch and
  entity types.
 */
template <class P> class tree_topology : public P, public data::data_client_t {

  /**
   * @ brief enum for the MPI asynchrnous communications
   */
  enum mpi_comm : int {
    LOCAL_REQUEST = 1,
    SOURCE_REQUEST = 2,
    MPI_DONE = 3,
    SOURCE_REPLY = 4,
    MPI_RANK_DONE = 5,
    FAILED_PROBE = 6,
    MPI_SHARE_EDGE_K = 7,
    MPI_SHARE_EDGE = 8
  };

public:
  using Policy = P; // Tree policy defined by the user

  static const size_t dimension = Policy::dimension; // Current dimension: 1,2,3
  using element_t = typename Policy::element_t; // Type of element either F or D
  using point_t = point__<element_t, dimension>;
  using range_t = std::array<point_t, 2>;
  using key_t = typename Policy::key_t;
  using branch_id_t = key_t;
  using branch_id_vector_t = std::vector<branch_id_t>;
  using branch_t = typename Policy::branch_t;
  using branch_vector_t = std::vector<branch_t *>;
  using entity_t = typename Policy::entity_t;
  using tree_entity_t = tree_entity<dimension, element_t, key_t, entity_t>;
  using entity_vector_t = std::vector<tree_entity_t *>;
  using apply_function = std::function<void(branch_t &)>;
  using entity_id_vector_t = std::vector<size_t>;
  using geometry_t = tree_geometry<element_t, dimension>;
  using entity_space_ptr_t = std::vector<tree_entity_t *>;

  // Hasher for the branch id used in the unordered_map data structure
  template <class KEY> struct branch_id_hasher__ {
    size_t operator()(const KEY &k) const noexcept {
      return k.value() & ((1 << 22) - 1);
    }
  };

  /*!
    Constuct a tree topology with unit coordinates, i.e. each coordinate
    dimension is in range [0, 1].
   */
  tree_topology() {
    branch_map_.emplace(branch_id_t::root(), branch_id_t::root());
    root_ = branch_map_.find(branch_id_t::root());
    assert(root_ != branch_map_.end());

    max_depth_ = 0;
    ghosts_entities_.resize(max_traversal);
    current_ghosts = 0;
  }

  /*!
    Construct a tree topology with specified ranges [end, start] for each
    dimension.
   */
  tree_topology(const point__<element_t, dimension> &start,
                const point__<element_t, dimension> &end) {
    branch_map_.emplace(branch_id_t::root(), branch_id_t::root());
    root_ = branch_map_.find(branch_id_t::root());
    assert(root_ != branch_map_.end());

    max_depth_ = 0;
    range_[0] = start;
    range_[1] = end;
    ghosts_entities_.resize(max_traversal);
    current_ghosts = 0;
  }

  /**
   * @brief Destroy the tree: empty the hash-table and destroy the entities
   * lists
   */
  ~tree_topology() {
    branch_map_.clear();
    tree_entities_.clear();
    ghosts_entities_.clear();
    current_ghosts = 0;
    shared_entities_.clear();
  }

  /**
   * Clean the tree topology but not the local bodies
   */
  void clean() {
    branch_map_.clear();
    tree_entities_.clear();
    for (int i = 0; i <= current_ghosts; ++i)
      ghosts_entities_[i].clear();
    current_ghosts = 0;
    shared_entities_.clear();
    branch_map_.emplace(branch_id_t::root(), branch_id_t::root());
    root_ = branch_map_.find(branch_id_t::root());
    assert(root_ != branch_map_.end());
    max_depth_ = 0;
  }

  /**
   * Reset the ghosts local information for the next tree traversal
   */
  void reset_ghosts(bool do_share_edge = true) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
      current_ghosts = 0;
      return;
    }
    // Same for the shared_edge entities
    // More complex because we need to find them and remove elements
    for (auto &s : shared_entities_) {
      // Find the parent and change the status
      key_t k = s.key();
      auto &b = find_parent(k);
      std::vector<size_t> remove;
      for (auto eid : b) {
        auto ent = get(eid);
        if (ent->owner() != rank) {
          ent->set_entity_ptr(nullptr);
          remove.push_back(eid);
        }
      }
      for (auto r : remove) {
        b.remove(r);
      }
      b.set_ghosts_local(false);
      b.set_requested(false);
      if (b.is_shared()) {
        b.set_locality(branch_t::LOCAL);
      }
    }

    for (int i = 0; i <= current_ghosts; ++i) {
      for (auto &g : ghosts_entities_[i]) {
        // Find the parent and change the status
        key_t k = g.key();
        auto &b = find_parent(k);
        assert(!b.is_local());
        for (auto eid : b) {
          auto ent = get(eid);
          assert(ent->owner() != rank);
          ent->set_entity_ptr(nullptr);
        }
        b.clear();
        b.set_ghosts_local(false);
        b.set_requested(false);
      }
    }

    // Find the first position of ghosts in tree_entities_
    // Empty the ghosts arrays
    size_t index_ghosts = tree_entities_.size();
    for (int i = 0; i <= current_ghosts; ++i) {
      index_ghosts -= ghosts_entities_[i].size();
      ghosts_entities_[i].clear();
    }
    current_ghosts = 0;
    index_ghosts -= shared_entities_.size();
    shared_entities_.clear();
    assert(!tree_entities_[index_ghosts].is_local());

    tree_entities_.erase(tree_entities_.begin() + index_ghosts,
                         tree_entities_.end());
    // Do not share edge in the case of keeping the tree
    if (do_share_edge) {
      share_edge();
    } else {
      remove_non_local();
    }
    cofm(root(), 0, false);
  }

  /**
   * \brief Share the edge particles to my direct neighbors
   * regarding the key ordering: 0 <-> 1 <-> 2 <-> 3 for 4 processes
   */
  /**
   * \brief Share the edge particles to my direct neighbors
   * regarding the key ordering: 0 <-> 1 <-> 2 <-> 3 for 4 processes
   */
  void share_edge() {
    // Communications 2 by 2
    // Send my highest key and my lowest key
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // First rank%2 == 0 send to rank%2 == 1
    std::array<key_t, 2> my_keys;
    std::array<key_t, 2> neighbor_keys;
    size_t byte_size = sizeof(std::array<key_t, 2>);

    std::vector<entity_t> received_ghosts;

    int partner = rank;
    for (int i = 0; i < 2; ++i) {
      if (rank % 2 == i) {
        partner = rank + 1;
        if (partner < size) {
          // Get my last key
          my_keys[0] = entities_.back().key();
          // Get the parent of this entity
          my_keys[1] = find_parent(my_keys[0]).key();
          MPI_Send(&(my_keys[0]), byte_size, MPI_BYTE, partner,
              MPI_SHARE_EDGE_K,MPI_COMM_WORLD);
          MPI_Recv(&(neighbor_keys[0]), byte_size, MPI_BYTE, partner,
              MPI_SHARE_EDGE_K,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      } else {
        partner = rank - 1;
        if (partner >= 0) {
          // Get my last key
          my_keys[0] = entities_.front().key();
          // Get the parent of this entity
          my_keys[1] = find_parent(my_keys[0]).key();
          MPI_Recv(&(neighbor_keys[0]), byte_size, MPI_BYTE, partner,
              MPI_SHARE_EDGE_K,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Send(&(my_keys[0]), byte_size, MPI_BYTE, partner,
              MPI_SHARE_EDGE_K,MPI_COMM_WORLD);
        }
      }
      // Entities to send to other rank
      std::vector<entity_t> edge_entities;
      if (partner >= 0 && partner < size) {
        // Compute the conflict
        int neighbor_parent_depth = neighbor_keys[1].depth();
        int my_parent_depth = my_keys[1].depth();
        // Handle the cases
        // The keys of the parents are the same
        if (neighbor_keys[1] == my_keys[1]) {
          // Send the bodies with this key's begining
          // If this rank is the small one, go from back of vector
          int max_entities = 1 << dimension;
          if (rank < partner) {
            auto cur = entities_.rbegin();
            auto end = entities_.rend();
            for (; cur != end && max_entities > 0; ++cur, --max_entities) {
              auto b = cur;
              key_t key = cur->key();
              key.truncate(neighbor_parent_depth);
              if (key == neighbor_keys[1]) {
                edge_entities.push_back(*b);
              }
            }
            // If this rank is the highest, go from the front of the vector
          } else {
            auto cur = entities_.begin();
            auto end = entities_.end();
            for (; cur != end && max_entities > 0; ++cur, --max_entities) {
              auto b = cur;
              key_t key = cur->key();
              key.truncate(neighbor_parent_depth);
              if (key == neighbor_keys[1]) {
                edge_entities.push_back(*b);
              }
            }
          }
        } else { // Send nothing if the parent are not the same, other branch
        }
        // If the rank has the highest branch
        // In this case just look if my body go in the neighbor branch
        if (neighbor_parent_depth < my_parent_depth) {
          // Is there a conflict: in this case compare the bodies
          key_t nb_key = neighbor_keys[0];
          nb_key.truncate(my_parent_depth);
          key_t my_key = my_keys[0];
          my_key.truncate(my_parent_depth);
          // I send the bodies that go in conflict with its body
          if (nb_key == my_key) {
            int max_entities = 1 << dimension;
            if (rank < partner) {
              auto cur = entities_.rbegin();
              auto end = entities_.rend();
              for (; cur != end && max_entities > 0; ++cur, --max_entities) {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(my_parent_depth);
                if (key == my_keys[1]) {
                  edge_entities.push_back(*b);
                }
              }
              // If this rank is the highest, go from the front of the vector
            } else {
              auto cur = entities_.begin();
              auto end = entities_.end();
              for (; cur != end && max_entities > 0; ++cur, --max_entities) {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(my_parent_depth);
                if (key == my_keys[1]) {
                  edge_entities.push_back(*b);
                }
              }
            }
          } else {
          }
        }
        // If my neighbor have the highest branch
        if (neighbor_parent_depth > my_parent_depth) {

          // Refine to the neighbor's depth
          while (my_parent_depth < neighbor_parent_depth) {
            refine_(find_parent(my_keys[0]));
            ++my_parent_depth;
          }

          key_t my_key = my_keys[0];
          my_key.truncate(neighbor_parent_depth);
          if (my_key == neighbor_keys[1]) {
            // Send my branch entities
            int max_entities = 1 << dimension;
            if (rank < partner) {
              auto cur = entities_.rbegin();
              auto end = entities_.rend();
              for (; cur != end && max_entities > 0; ++cur, --max_entities) {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(neighbor_parent_depth);
                if (key == neighbor_keys[1]) {
                  edge_entities.push_back(*b);
                }
              }
              // If this rank is the highest, go from the front of the vector
            } else {
              auto cur = entities_.begin();
              auto end = entities_.end();
              for (; cur != end && max_entities > 0; ++cur, --max_entities) {
                auto b = cur;
                key_t key = cur->key();
                key.truncate(neighbor_parent_depth);
                if (key == neighbor_keys[1]) {
                  edge_entities.push_back(*b);
                }
              }
            }
          } else {
            // Nothing to share but be prepare to receive
            // the other bodies at the right depth
            int nrefine = neighbor_parent_depth - my_parent_depth;
            int depth = my_parent_depth;
            for (int r = 0; r < nrefine; ++r) {
              // Find my branch
              key_t branch_key = my_keys[0];
              branch_key.truncate(my_parent_depth);
              auto &b = branch_map_.find(branch_key)->second;
              refine_(b);
              ++my_parent_depth;
            }
          }
        }
      }
      // 3. Send them
      if (rank % 2 == i) {
        partner = rank + 1;
        if (partner < size) {
          MPI_Send(&(edge_entities[0]), sizeof(entity_t) * edge_entities.size(),
                   MPI_BYTE, partner, MPI_SHARE_EDGE, MPI_COMM_WORLD);
          MPI_Status status;
          MPI_Probe(partner, MPI_SHARE_EDGE, MPI_COMM_WORLD, &status);
          // Get the count
          int nrecv = 0;
          MPI_Get_count(&status, MPI_BYTE, &nrecv);
          int offset = received_ghosts.size();
          received_ghosts.resize(offset + nrecv / sizeof(entity_t));
          MPI_Recv(&(received_ghosts[offset]), nrecv, MPI_BYTE, partner, MPI_SHARE_EDGE,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      } else {
        partner = rank - 1;
        if (partner >= 0) {
          MPI_Status status;
          MPI_Probe(partner, MPI_SHARE_EDGE, MPI_COMM_WORLD, &status);
          // Get the count
          int nrecv = 0;
          MPI_Get_count(&status, MPI_BYTE, &nrecv);
          int offset = received_ghosts.size();
          received_ghosts.resize(offset + nrecv / sizeof(entity_t));
          MPI_Recv(&(received_ghosts[offset]), nrecv, MPI_BYTE, partner, MPI_SHARE_EDGE,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Send(&(edge_entities[0]), sizeof(entity_t) * edge_entities.size(),
                   MPI_BYTE, partner, MPI_SHARE_EDGE, MPI_COMM_WORLD);
        }
      } // if
    } // for

    for (auto g : received_ghosts) {
      shared_entities_.push_back(g);
    }
    for (auto &g : shared_entities_) {
      auto *bi = &(g);
      assert(bi->mass() != 0.);
      assert(bi->key().value() != 0);
      if (bi->key() > entities().back().key()) {
        partner = rank + 1;
      } else {
        partner = rank - 1;
      }
      auto id = make_entity(bi->key(), bi->coordinates(), nullptr, partner,
                            bi->mass(), bi->id(), bi->radius());
      insert(id);
      get(id)->set_entity_ptr(bi);
      // Set the ghosts local in this case
    }
    for (auto &g : shared_entities_) {
      find_parent(g.key()).set_ghosts_local(true);
    }
  }

  /**
   * \brief Change the range of the tree topology
   */
  void set_range(const range_t &range) { range_ = range; }

  /**
   * @brief Get the range
   */
  const std::array<point__<element_t, dimension>, 2> &range() { return range_; }

  /**
   * @brief Get the ci-th child of the given branch.
   */
  branch_t *child(branch_t *b, size_t ci) {
    // Use the hash table
    branch_id_t bid = b->key();         // Branch id
    bid.push(ci);                       // Add child number
    auto child = branch_map_.find(bid); // Search for the child
    // If it does not exists, return nullptr
    if (child == branch_map_.end())
      return nullptr;
    return &child->second;
  }

  /**
   * @brief Return reference to the vector of tree_entities
   */
  std::vector<tree_entity_t> &tree_entities() { return tree_entities_; }

  /**
   * @ brief Return a reference to the vector of the entities
   */
  std::vector<entity_t> &entities() { return entities_; }

  /**
   * @brief Return an entity by its id
   */
  template <typename E> entity_t &entity(E e) {
    return entities_[static_cast<int>(e)];
  }

  /**
   * @brief Return a vector with all the local sub entities
   *
   * @param start The branch in which the search occur
   * @param search_list The found entities vector
   */
  void get_sub_entities(branch_t *start, std::vector<entity_t> &search_list) {
    std::stack<branch_t *> stk;
    stk.push(start);

    while (!stk.empty()) {
      branch_t *c = stk.top();
      stk.pop();
      if (c->is_leaf()) {
        for (auto id : *c) {
          auto child = this->get(id);
          if (child->entity_ptr() != nullptr) {
            search_list.push_back(*(child->entity_ptr()));
          }
        }
      } else {
        for (int i = 0; i < (1 << dimension); ++i) {
          if (!c->as_child(i))
            continue;
          auto next = child(c, i);
          assert(next != nullptr);
          stk.push(next);
        }
      }
    }
  }

  /**
   * @brief Find all the center of mass of the tree up to the
   * maximum sub particles criterion.
   *
   * @param b The starting branch for the search, usually root
   * @param criterion The maximum of subparticles for the COMs
   * @param search_list The extracted COMs
   */
  void find_sub_cells(branch_t *b, const uint64_t &criterion,
                      std::vector<branch_t *> &search_list) {
    std::stack<branch_t *> stk;
    stk.push(b);

    while (!stk.empty()) {
      branch_t *c = stk.top();
      stk.pop();
      if (c->is_leaf() && c->is_local()) {
        search_list.push_back(c);
      } else {
        if (c->is_local() && c->sub_entities() <= criterion) {
          search_list.push_back(c);
        } else {
          for (int i = (1 << dimension) - 1; i >= 0; --i) {
            if (!c->as_child(i))
              continue;
            auto next = child(c, i);
            stk.push(next);
          }
        }
      }
    }
  }

  /**
` * @brief Apply a function ef to the sub_cells using asynchronous comms.
  * @param [in] b The starting branch for the traversal
  * @param [in] ncritical The number of subparticles for the branches
  * @param [in] b The starting branch for the traversal
  * @param [in] b The starting branch for the traversal
  * @return <return_description>
  * @details <details>
  */
  template <typename EF, typename... ARGS>
  void traversal_sph(branch_t *b, EF &&ef, ARGS &&... args) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    entities_w_ = entities_;

    std::vector<branch_t *> working_branches;
    find_sub_cells(b, ncritical, working_branches);

    // Remaining branches in case of non locality
    std::vector<branch_t *> remaining_branches;

    // Start communication thread
    if (size != 1) {
      std::thread handler(&tree_topology::handle_requests, this);
      // Start tree traversal
      traverse_sph(working_branches, remaining_branches, false, ef,
                   std::forward<ARGS>(args)...);
      MPI_Send(NULL, 0, MPI_INT, rank, MPI_DONE, MPI_COMM_WORLD);
      // Wait for communication thread
      handler.join();
    } else {
      traverse_sph(working_branches, remaining_branches, false, ef,
                   std::forward<ARGS>(args)...);
    }

#ifdef DEBUG
    int flag = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
               &status);
    if (flag != 0) {
      std::cerr <<rank<< " TAG: " << status.MPI_TAG << " SOURCE: "<<
        status.MPI_SOURCE<< std::endl;
    };
    assert(flag == 0);
#endif

    // Add the eventual ghosts in the tree for remaining branches
    for (size_t i = 0; i < ghosts_entities_[current_ghosts].size(); ++i) {
      entity_t &g = ghosts_entities_[current_ghosts][i];
      auto id = make_entity(g.key(), g.coordinates(), nullptr, g.owner(),
                            g.mass(), g.id(), g.radius());
      // Assert the parent exists and is non local
      insert(id);
      auto nbi = get(id);
      nbi->set_entity_ptr(&g);
      find_parent(g.key()).set_ghosts_local(true);
    }

    // Prepare for the eventual next tree traversal, use other ghosts vector
    ++current_ghosts;
    assert(current_ghosts < max_traversal);

    // Recompute the COFM because new entities are present in the tree
    // Vector useless because in this case no ghosts can be found
    if (remaining_branches.size() > 0) {
      std::vector<branch_t *> ignore;
      traverse_sph(remaining_branches, ignore, true, ef,
                   std::forward<ARGS>(args)...);
    }
    // Copy back the results
    entities_ = entities_w_;
    // Synchronize the threads to be sure they dont start to share edges
    // Critical
    MPI_Barrier(MPI_COMM_WORLD);

  } // apply_sub_cells

  /**
   * @brief Perform a tree traversal in parallel using omp threads for
   * each working branch. If a branch is not local, make a request to the
   * thread handler.
   * @param [in] do_square True is varible smoothing length false otherwize
   * @param [in] work_branch The set of branch on which to compute the inter.
   * @param [in] non_local_branches The branches identified as non local
   * during this traversal, they will be computed after this traversal when the
   * distant particles will be gathered
   * @param [in] assert_local In the case of local run assert that no distant
   * particles are reached.
   * @param [in] ef The function to apply to the particles in the work_branch
   * @param [in] args Arguments of the function ef
   * @return void
   * @details
   */
  template <typename EF, typename... ARGS>
  void traverse_sph(std::vector<branch_t *> &working_branches,
                    std::vector<branch_t *> &non_local_branches,
                    const bool assert_local, EF &&ef, ARGS &&... args) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nelem = working_branches.size();

    size_t max_send = 20;

#pragma omp parallel
    {
      std::vector<branch_t *> omp_non_local;
      std::vector<key_t> send;
#pragma omp for schedule(static)
      for (int i = 0; i < nelem; ++i) {
        std::vector<branch_t *> inter_list;
        std::vector<branch_t *> requests_branches;

        // Compute the interaction list for this branch
        if (interactions_branches(working_branches[i], inter_list,
                                  requests_branches)) {
          std::vector<std::vector<entity_t *>> neighbors;
          neighbors.clear();
          neighbors.resize(working_branches[i]->sub_entities());
          // Sub traversal to apply to the particles
          interactions_particles(working_branches[i], inter_list, neighbors);
          // Perform the computation in the same time for all the threads
          int index = 0;
          // for(auto j: *(working_branches[i]))
          for (int j = working_branches[i]->begin_tree_entities();
               j <= working_branches[i]->end_tree_entities(); ++j) {
            if (tree_entities_[j].is_local())
              ef(entities_w_[j], neighbors[index], std::forward<ARGS>(args)...);
            ++index;
          }
        } else {
          assert(!assert_local);
          omp_non_local.push_back(working_branches[i]);
#pragma omp critical
          {
            // Send branch key to request handler
            for (auto b : requests_branches) {
              if (!b->requested()) {
                send.push_back(b->key());
                b->set_requested(true);
              }
            }
          } // omp critical
          if (send.size() >= max_send) {
            MPI_Send(&(send[0]), send.size() * sizeof(key_t), MPI_BYTE, rank,
                     LOCAL_REQUEST, MPI_COMM_WORLD);
            // Reset send
            send.clear();
          }
        } // if else
      }   // for
      // Finish the send
      if (send.size() > 0) {
        MPI_Send(&(send[0]), send.size() * sizeof(key_t), MPI_BYTE, rank,
                 LOCAL_REQUEST, MPI_COMM_WORLD);
        // Reset send
        send.clear();
      }
#pragma omp critical
      non_local_branches.insert(non_local_branches.begin(),
                                omp_non_local.begin(), omp_non_local.end());
    } // omp parallel
  }   // traverse_sph

  /**
   * @brief Compute the interaction list for all the sub-particles in b
   * @param [in] b Branch on which to propagate the force
   * @param [in] inter_list The interaction list of the sub-particles of b
   * @param [in] ef The function to apply on the particles
   * @param [in] args The arguments of the function to apply
   * @return void
   * @details This function apply a tree traversal because the branch b can
   * be different than a leaf.
   */
  bool interactions_branches(branch_t *work_branch,
                             std::vector<branch_t *> &inter_list,
                             std::vector<branch_t *> &non_local) {
    // Queues for the branches
    std::vector<branch_t *> queue;
    std::vector<branch_t *> new_queue;

    bool missing_branch = false;

    queue.push_back(root());
    while (!queue.empty()) {
      new_queue.clear();
      // Add the next level in the queue
      for (int i = 0; i < queue.size(); ++i) {
        branch_t *c = queue[i];
        assert(!c->is_leaf());
        for (int d = 0; d < (1 << dimension); ++d) {
          if (!c->as_child(d))
            continue;
          new_queue.push_back(child(c, d));
        }
      }
      const int queue_size = new_queue.size();
      queue.clear();
      for (int i = 0; i < queue_size; ++i) {
        branch_t *b = new_queue[i];
        if (geometry_t::intersects_box_box(b->bmin(), b->bmax(),
                                           work_branch->bmin(),
                                           work_branch->bmax())) {
          if (b->is_leaf()) {
            if (b->is_local() || b->ghosts_local()) {
              inter_list.push_back(b);
            } else {
              missing_branch = true;
              if (!b->requested())
                non_local.push_back(b);
            }
          } else {
            queue.push_back(new_queue[i]);
          }
        }
      }
    }
    return !missing_branch;
  }

  void interactions_particles(branch_t *working_branch,
                              const std::vector<branch_t *> &inter_list,
                              std::vector<std::vector<entity_t *>> &neighbors) {
    std::vector<point_t> inter_coordinates;
    std::vector<element_t> inter_radius;
    std::vector<entity_t *> inter_entities;
    for (int j = 0; j < inter_list.size(); ++j) {
      for (auto k : *(inter_list[j])) {
        inter_coordinates.push_back(tree_entities_[k].coordinates());
        inter_radius.push_back(tree_entities_[k].radius());
        inter_entities.push_back(tree_entities_[k].entity_ptr());
        assert(inter_entities.back() != nullptr);
      }
    }
    const int nb_entities = inter_coordinates.size();
    int index = 0;
    for (int i = working_branch->begin_tree_entities();
         i <= working_branch->end_tree_entities(); ++i) {
      point_t coordinates = tree_entities_[i].coordinates();
      element_t radius = tree_entities_[i].radius();
      size_t total = 0;
      std::vector<size_t> accepted(nb_entities, 0);
      for (int j = 0; j < nb_entities; ++j) {
        accepted[j] += geometry_t::within_square(
            inter_coordinates[j], coordinates, inter_radius[j], radius);
        total += accepted[j];
      }
      neighbors[index].resize(total);
      int index_add = 0;
      for (int j = 0; j < nb_entities; ++j) {
        if (accepted[j])
          neighbors[index][index_add++] = inter_entities[j];
      }
      ++index;
    }
  }

  template <typename FC, typename DFCDR, typename DFCDRDR, typename C2P>
  void traversal_fmm(branch_t *b, double maxmasscell, const double MAC,
                     FC &&f_fc, DFCDR &&f_dfcdr, DFCDRDR &&f_dfcdrdr,
                     C2P &&f_c2p) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    clog_one(trace) << "FMM : maxmasscell: " << maxmasscell << " MAC: " << MAC
                    << std::endl;

    entities_w_ = entities_;

    std::vector<branch_t *> work_branch;
    std::vector<branch_t *> remaining_branches;
    find_sub_cells(b, 32, work_branch);

    if (size != 1) {
      std::thread handler(&tree_topology::handle_requests, this);
      // Start tree traversal
      traverse_fmm(work_branch, remaining_branches, MAC, f_fc, f_dfcdr,
                    f_dfcdrdr, f_c2p);
      MPI_Send(NULL, 0, MPI_INT, rank, MPI_DONE, MPI_COMM_WORLD);
      // Wait for communication thread
      handler.join();
    } else {
      traverse_fmm(work_branch, remaining_branches, MAC, f_fc, f_dfcdr,
                    f_dfcdrdr, f_c2p);
    }

    // Check if no message remainig
#ifdef DEBUG
    int flag = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
               &status);
    if(flag != 0){
      std::cerr<<rank<<" TAG:"<<status.MPI_TAG<<" SOURCE: "<<
        status.MPI_SOURCE<<std::endl;
    }
    assert(flag == 0);
#endif
    for (size_t i = 0; i < ghosts_entities_[current_ghosts].size(); ++i) {
      entity_t &g = ghosts_entities_[current_ghosts][i];
      auto id = make_entity(g.key(), g.coordinates(), nullptr, g.owner(),
                            g.mass(), g.id(), g.radius());
      insert(id);
      auto nbi = get(id);
      nbi->set_entity_ptr(&g);
      // Set the parent to local for the search
      find_parent(g.key()).set_ghosts_local(true);
    }
    ++current_ghosts;
    assert(current_ghosts < max_traversal);

    // compute the particle to particle interaction on each sub-branches
    if (size != 1) {
      std::vector<branch_t *> ignore;
      traverse_fmm(remaining_branches, ignore, MAC, f_fc, f_dfcdr, f_dfcdrdr,
                    f_c2p);
    } else {
      assert(remaining_branches.size() == 0);
    }
    entities_ = entities_w_;


    // Synchronize the threads to be sure they dont start to share edges
    // Critical
    MPI_Barrier(MPI_COMM_WORLD);

  } // apply_sub_cells

  /**
   * @brief Similar to traverse_sph but using the MAC criterion and applying
   * the C2C and C2P computations
   */
  template <typename FC, typename DFCDR, typename DFCDRDR, typename C2P>
  void traverse_fmm(std::vector<branch_t *> &work_branch,
                     std::vector<branch_t *> &remaining_branches,
                     const double MAC, FC &&f_fc, DFCDR &&f_dfcdr,
                     DFCDRDR &&f_dfcdrdr, C2P &&f_c2p) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t total_done = 0;

    size_t nelem = work_branch.size();

#pragma omp parallel for
    for (size_t i = 0; i < nelem; ++i) {
      std::vector<branch_t *> inter_list;
      std::vector<branch_t *> requests_branches;

      total_done += work_branch[i]->sub_entities();
      // Sub traversal and compute the MAC

      if (traversal_c2c_c2p_p2p(work_branch[i], requests_branches, MAC, f_fc,
                                f_dfcdr, f_dfcdrdr, f_c2p)) {
        if(size == 1 ) assert(false);
        std::vector<key_t> send;
#pragma omp critical
        {
          remaining_branches.push_back(work_branch[i]);
          // Send branch key to request handler
          for (auto b : requests_branches) {
            assert(b->owner() < size && b->owner() >= 0);
            assert(b->owner() != rank);
            if (!b->requested()) {
              send.push_back(b->key());
              b->set_requested(true);
            }
          }
        }
        MPI_Request request;
        MPI_Send(&(send[0]), send.size() * sizeof(key_t), MPI_BYTE, rank,
                 LOCAL_REQUEST, MPI_COMM_WORLD);
      }
    }
  }

  /**
   * @brief Center of Mass (COM) to COM interaction and the COM to the particles
   * in this COM.
   */
  template <typename FC, typename DFCDR, typename DFCDRDR, typename F_C2P>
  bool traversal_c2c_c2p_p2p(branch_t *b, std::vector<branch_t *> &non_local,
                             const double MAC, FC &&f_fc, DFCDR &&f_dfcdr,
                             DFCDRDR &&f_dfcdrdr, F_C2P &&f_c2p) {
    point_t fc = 0.;
    double dfcdr[9] = {0.};
    double dfcdrdr[27] = {0.};
    point_t coordinates = b->coordinates();

    std::vector<branch_t *> queue;
    std::vector<branch_t *> new_queue;
    std::vector<point_t> c2c_coordinates;
    std::vector<double> c2c_mass;

    // -------- 1. C2C interactions, keep track of the leaves non interacted
    std::vector<branch_t *> interactions_leaves;
    queue.push_back(root());
    while (!queue.empty()) {
      new_queue.clear();
      const int queue_size = queue.size();
      for (int i = 0; i < queue_size; ++i) {
        branch_t *q = queue[i];
        for (int d = 0; d < (1 << dimension); ++d) {
          if (!q->as_child(d))
            continue;
          branch_t *c = child(q, d);
          if (geometry_t::box_MAC(coordinates, c->coordinates(), c->bmin(),
                                  c->bmax(), MAC)) {
            c2c_coordinates.push_back(c->coordinates());
            c2c_mass.push_back(c->mass());
          } else if (!c->is_leaf()) {
            new_queue.push_back(c);
          } else if (!c->is_local() && !c->ghosts_local()) {
            non_local.push_back(c);
          } else {
            interactions_leaves.push_back(c);
          }
        }
      }
      queue.clear();
      queue = new_queue;
    }

    // Compute the C2C
    for (int i = 0; i < c2c_coordinates.size(); ++i) {
      // Compute the matrices
      f_fc(fc, coordinates, c2c_coordinates[i], c2c_mass[i]);
      f_dfcdr(dfcdr, coordinates, c2c_coordinates[i], c2c_mass[i]);
      f_dfcdrdr(dfcdrdr, coordinates, c2c_coordinates[i], c2c_mass[i]);
    }

    // If all the sub particles are present
    if (non_local.size() == 0) {
      // Propagate this information to the sub-particles for C2P
      for (int i = b->begin_tree_entities(); i <= b->end_tree_entities(); ++i) {
        if (tree_entities_[i].is_local()) {
          f_c2p(fc, dfcdr, dfcdrdr, b->coordinates(), &(entities_w_[i]));
        }
      }
      // Apply the P2P for the local particles
      for (int l = 0; l < interactions_leaves.size(); ++l) {
        for (auto k : *(interactions_leaves[l])) {
          for (int i = b->begin_tree_entities(); i <= b->end_tree_entities();
               ++i) {
            if (tree_entities_[k].id() != tree_entities_[i].id()) {
              // N square computation
              point_t fc;
              entities_w_[i].setAcceleration(
                  entities_w_[i].getAcceleration() +
                  f_fc(fc, entities_w_[i].coordinates(),
                       tree_entities_[k].coordinates(),
                       tree_entities_[k].mass()));
            }
          }
        }
      }
    }
    return non_local.size() > 0;
  }

  /**
   *
   */
  void handle_requests() {
    const int max_size = 500;
    const int max_requests = 1000;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    bool done_traversal = false;
    bool done = false;

    // Maintain a request array for all neighbors
    std::vector<std::vector<std::vector<key_t>>> requests(size);
    std::vector<std::vector<std::vector<entity_t>>> reply(size);
    std::vector<int> current_requests(size, 0);
    std::vector<int> current_reply(size, 0);
    for (int i = 0; i < size; ++i) {
      requests[i].resize(max_requests);
      reply[i].resize(max_requests);
    }

    std::vector<MPI_Request> mpi_requests(max_requests);
    int current_mpi_requests = 0;
    std::vector<MPI_Request> mpi_replies(max_requests);
    int current_mpi_replies = 0;
    int request_counter = 0;
    std::vector<bool> rank_done(size, false);

    while (!done) {
      MPI_Status status;
      // Wait on probe
      int source, tag, nrecv;
      if (!done_traversal) {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        source = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        nrecv = 0;
        MPI_Get_count(&status, MPI_BYTE, &nrecv);
      } else {
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
          source = status.MPI_SOURCE;
          tag = status.MPI_TAG;
          nrecv = 0;
          MPI_Get_count(&status, MPI_BYTE, &nrecv);
        } else {
          tag = FAILED_PROBE;
        }
      }

      switch (tag) {
      case MPI_RANK_DONE: {
        MPI_Recv(NULL, 0, MPI_INT, source, MPI_RANK_DONE, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        rank_done[source] = true;
      } break;
      // ------------------------------------------------------------------ //
      //     Another rank replied to my entities request                    //
      // ------------------------------------------------------------------ //
      case SOURCE_REPLY: {
        assert(nrecv != 0);
        int current = ghosts_entities_[current_ghosts].size();
        ghosts_entities_[current_ghosts].resize(current +
                                                nrecv / sizeof(entity_t));
        MPI_Recv(&(ghosts_entities_[current_ghosts][current]), nrecv, MPI_BYTE,
                 source, SOURCE_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        --request_counter;
      } break;
      // ------------------------------------------------------------------ //
      //    Another rank request information for entities                   //
      // ------------------------------------------------------------------ //
      case SOURCE_REQUEST: {
        std::vector<key_t> received(nrecv / sizeof(key_t));
        MPI_Recv(&(received[0]), nrecv, MPI_BYTE, source, SOURCE_REQUEST,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (auto k : received) {
          auto branch = branch_map_.find(k);
          assert(branch != branch_map_.end());
          get_sub_entities(&(branch->second),
                           reply[source][current_reply[source]]);
        }
        int ncount = reply[source][current_reply[source]].size();
        MPI_Isend(&(reply[source][current_reply[source]++][0]),
                  ncount * sizeof(entity_t), MPI_BYTE, source, SOURCE_REPLY,
                  MPI_COMM_WORLD, &(mpi_replies[current_mpi_replies++]));
        if (current_mpi_replies > max_requests) {
          clog_one(error) << rank << ": Exceeding number of replies requests"
                          << std::endl;
        }
      } break;
      // ------------------------------------------------------------------ //
      //        A local thread requested a distant particles                //
      // ------------------------------------------------------------------ //
      case LOCAL_REQUEST: {
        assert(done_traversal == false);
        std::vector<key_t> keys(nrecv / sizeof(key_t));
        MPI_Recv(&(keys[0]), nrecv, MPI_BYTE, rank, LOCAL_REQUEST,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (auto k : keys) {
          auto itr = branch_map_.find(k);
          assert(itr != branch_map_.end());
          branch_t *branch = &(itr->second);
          assert(branch->requested());
          int owner = branch->owner();
          assert(owner != rank);
          // Add this request to vector
          requests[owner][current_requests[owner]].push_back(k);
          if (requests[owner][current_requests[owner]].size() >= max_size) {
            MPI_Isend(&(requests[owner][current_requests[owner]++][0]),
                      max_size * sizeof(key_t), MPI_BYTE, owner, SOURCE_REQUEST,
                      MPI_COMM_WORLD, &(mpi_requests[current_mpi_requests++]));
            if (current_mpi_requests > max_requests) {
              clog_one(error)
                  << rank << ": Exceeding number of requests requests"
                  << std::endl;
            }
            ++request_counter;
          }
        }
      } break;
      // ------------------------------------------------------------------ //
      //        The OpenMP threads are done, send the last requests         //
      // ------------------------------------------------------------------ //
      case MPI_DONE: {
        // First time, send remaining requests
        MPI_Recv(NULL, 0, MPI_INT, source, MPI_DONE, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        for (int i = 0; i < size; ++i) {
          if (i == rank)
            continue;
          if (requests[i][current_requests[i]].size() > 0) {
            int nsend = requests[i][current_requests[i]].size();
            MPI_Isend(&(requests[i][current_requests[i]][0]),
                      nsend * sizeof(key_t), MPI_BYTE, i, SOURCE_REQUEST,
                      MPI_COMM_WORLD, &(mpi_requests[current_mpi_requests++]));
            if (current_mpi_requests > max_requests) {
              clog_one(error)
                  << rank << ": Exceeding number of requests requests"
                  << std::endl;
            }
            ++request_counter;
          }
        }
        done_traversal = true;
      } break;

      case FAILED_PROBE:
        break;

      default: {
        std::cerr<<rank<<" TAG: "<<tag<<" FROM: "<<source<<std::endl;
        assert(false);
      } break;
      }

      if (done_traversal) {
        // Check if all requests have been answered
        done = request_counter == 0;
        bool done_requests = true;
        int flag;
        for (size_t i = 0; i < current_mpi_requests; ++i) {
          MPI_Test(&(mpi_requests[i]), &flag, MPI_STATUS_IGNORE);
          done_requests = done_requests && flag;
        }
        bool done_replies = true;
        for (size_t i = 0; i < current_mpi_replies; ++i) {
          MPI_Test(&(mpi_replies[i]), &flag, MPI_STATUS_IGNORE);
          done_replies = done_replies && flag;
        }
        done = done && done_replies && done_requests;
        // Just do this step once to send the last requests
        if (done && (!rank_done[rank])) {
          rank_done[rank] = true;
          // Send rank done
          for (size_t i = 0; i < size; ++i) {
            if (i == rank)
              continue;
            MPI_Isend(NULL, 0, MPI_INT, i, MPI_RANK_DONE, MPI_COMM_WORLD,
                      &(mpi_requests[current_mpi_requests++]));
            if (current_mpi_requests > max_requests) {
              clog_one(error)
                  << rank << ": Exceeding number of requests requests"
                  << std::endl;
            }
          }
        }
        for (size_t i = 0; i < size; ++i) {
          done = done && rank_done[i];
        }
      }
    }
  }

  void find_level(branch_t *start, const int &level,
                  std::vector<branch_t *> &find) {
    std::stack<branch_t *> stk;
    stk.push(root());
    while (!stk.empty()) {
      branch_t *c = stk.top();
      int cur_level = c->key().depth();
      stk.pop();
      if (c->is_leaf() && c->is_local()) {
        find.push_back(c);
      } else {
        if (c->is_local() && cur_level == level) {
          find.push_back(c);
        } else {
          for (int i = (1 << dimension) - 1; i >= 0; --i) {
            if (!c->as_child(i))
              continue;
            auto next = child(c, i);
            stk.push(next);
          }
        }
      }
    }
  }

  void cofm(branch_t *start, element_t epsilon = 0, bool local = false) {
    // Find the sub particles on which we want to work
    std::vector<branch_t *> working_branches;
    std::stack<branch_t *> stk_remaining;
    // in 3d: 8^3 branches maximum (512)
    int level = 5;
    std::stack<branch_t *> stk;
    stk.push(root());
    while (!stk.empty()) {
      branch_t *c = stk.top();
      int cur_level = c->key().depth();
      stk.pop();
      if (c->is_leaf() && c->is_local()) {
        working_branches.push_back(c);
      } else {
        if (c->is_local() && cur_level == level) {
          working_branches.push_back(c);
        } else {
          stk_remaining.push(c);
          for (int i = (1 << dimension) - 1; i >= 0; --i) {
            if (!c->as_child(i))
              continue;
            auto next = child(c, i);
            stk.push(next);
          }
        }
      }
    }

    // Work in parallel on the sub branches
    const int nwork = working_branches.size();

#pragma omp parallel for
    for (int b = 0; b < nwork; ++b) {
      // Find the leave in order in these sub branches
      std::stack<branch_t *> stk1;
      std::stack<branch_t *> stk2;
      stk1.push(working_branches[b]);
      while (!stk1.empty()) {
        branch_t *cur = stk1.top();
        stk1.pop();
        stk2.push(cur);
        // Push children to stk1
        if (!cur->is_leaf()) {
          for (int i = 0; i < (1 << dimension); ++i) {
            if (!cur->as_child(i))
              continue;
            branch_t *next = child(cur, i);
            stk1.push(next);
          }
        }
      }
      // Finish the highest part of the tree in serial
      while (!stk2.empty()) {
        branch_t *cur = stk2.top();
        stk2.pop();
        update_COM(cur, epsilon, local);
      }
    }
    // Finish the high part of the tree on one thread
    while (!stk_remaining.empty()) {
      branch_t *cur = stk_remaining.top();
      stk_remaining.pop();
      update_COM(cur, epsilon, local);
    }
  }

  // Functions for the tree traversal
  void update_COM(branch_t *b, element_t epsilon = element_t(0),
                  bool local_only = false) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    element_t mass = 0;
    point_t bmax{}, bmin{};
    // element_t radius = 0.;
    point_t coordinates{};
    uint64_t nchildren = 0;
    int owner = b->owner();
    size_t begin_te = tree_entities_.size();
    size_t end_te = 0;
    for (size_t d = 0; d < dimension; ++d) {
      bmax[d] = -DBL_MAX;
      bmin[d] = DBL_MAX;
    }
    bool full_nonlocal = true, full_local = true;
    if (b->is_leaf()) {
      // For local branches, compute the radius
      if (b->is_local()) {
        int start = -1;
        int end = -1;
        for (auto child : *b) {
          auto ent = get(child);
          if (ent->is_local()) {
            if (start == -1)
              start = child;
            end = child;
          }
          owner = ent->owner();
          if (local_only && !ent->is_local()) {
            continue;
          }
          if (ent->owner() == rank)
            full_nonlocal = false;
          else
            full_local = false;
          ++nchildren;
          element_t childmass = ent->mass();
          for (size_t d = 0; d < dimension; ++d) {
            bmax[d] = std::max(bmax[d], ent->coordinates()[d] + epsilon +
                                            ent->radius() / 2.);
            bmin[d] = std::min(bmin[d], ent->coordinates()[d] - epsilon -
                                            ent->radius() / 2.);
          }
          coordinates += childmass * ent->coordinates();
          mass += childmass;
        }
        if (mass > element_t(0))
          coordinates /= mass;
        begin_te = start;
        end_te = end;
      } else {
        // For non local particles use existing value from remote
        coordinates = b->coordinates();
        bmin = b->bmin();
        bmax = b->bmax();
        mass = b->mass();
        nchildren = b->sub_entities();
      }
      // Locality for leaves
      if (full_nonlocal && !full_local)
        b->set_owner(owner);
      if (b->owner() == rank && full_local)
        b->set_locality(branch_t::LOCAL);
      else if (b->owner() == rank && !full_local)
        b->set_locality(branch_t::SHARED);
      else
        b->set_locality(branch_t::NONLOCAL);

    } else {
      bool local = false;
      bool nonlocal = false;
      for (int i = 0; i < (1 << dimension); ++i) {
        auto branch = child(b, i);
        if (branch == nullptr)
          continue;
        nchildren += branch->sub_entities();
        mass += branch->mass();
        if (branch->locality() == branch_t::LOCAL)
          local = true;
        if (branch->locality() == branch_t::NONLOCAL)
          nonlocal = true;
        if (branch->locality() == branch_t::SHARED)
          local = nonlocal = true;
        if (branch->mass() > 0) {
          for (size_t d = 0; d < dimension; ++d) {
            bmax[d] = std::max(bmax[d], branch->bmax()[d]);
            bmin[d] = std::min(bmin[d], branch->bmin()[d]);
          }
        }
        coordinates += branch->mass() * branch->coordinates();

        begin_te = std::min(begin_te, branch->begin_tree_entities());
        end_te = std::max(end_te, branch->end_tree_entities());
      }
      if (mass > element_t(0))
        coordinates /= mass;
      if (local && nonlocal)
        b->set_locality(branch_t::SHARED);
      if (local && !nonlocal)
        b->set_locality(branch_t::LOCAL);
      if (!local && nonlocal)
        b->set_locality(branch_t::NONLOCAL);
    }
    b->set_sub_entities(nchildren);
    b->set_coordinates(coordinates);
    b->set_mass(mass);
    b->set_bmin(bmin);
    b->set_bmax(bmax);
    assert(nchildren != 0);
    b->set_begin_tree_entities(begin_te);
    b->set_end_tree_entities(end_te);
  }

  branch_t *find_branch(const key_t &key) {
    auto b = branch_map_.find(key);
    // assert(b != branch_map_.end());
    if (b == branch_map_.end())
      return nullptr;
    return &(b->second);
  }

  void find_children(const key_t &key, std::vector<branch_t *> &children) {
    auto b = &(branch_map_.find(key)->second);
    for (int d = 0; d < (1 << dimension); ++d) {
      if (!b->as_child(d))
        continue;
      children.push_back(child(b, d));
    }
  }

  /*!
    Return an index space containing all entities within the specified
    spheroid.
   */
  template <typename EF>
  entity_space_ptr_t find_in_radius(const point_t &center, element_t radius,
                                    EF &&ef) {
    entity_space_ptr_t ents;

    // ITERATIVE VERSION
    std::stack<branch_t *> stk;
    stk.push(root());

    while (!stk.empty()) {
      branch_t *b = stk.top();
      stk.pop();
      if (b->is_leaf()) {
        for (auto id : *b) {
          auto child = &(tree_entities_[id]);
          // Check if in radius
          if (ef(center, child->coordinates(), radius, child->radius())) {
            ents.push_back(child);
          }
        }
      } else {
        for (int i = 0; i < (1 << dimension); ++i) {
          if (!b->as_child(i))
            continue;
          auto branch = child(b, i);
          assert(branch != nullptr);
          if (geometry_t::intersects_sphere_box(branch->bmin(), branch->bmax(),
                                                center, radius)) {
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
  template <typename EF>
  entity_space_ptr_t find_in_box(const point_t &min, const point_t &max,
                                 EF &&ef) {
    entity_space_ptr_t ents;

    // ITERATIVE VERSION
    std::stack<branch_t *> stk;
    stk.push(root());

    while (!stk.empty()) {
      branch_t *b = stk.top();
      stk.pop();
      if (b->is_leaf()) {
        for (auto id : *b) {
          auto child = &(tree_entities_[id]);
          // Check if in box
          if (ef(min, max, child->coordinates(), child->radius())) {
            ents.push_back(child);
          }
        }
      } else {
        for (int i = 0; i < (1 << dimension); ++i) {
          if (!b->as_child(i))
            continue;
          auto branch = child(b, i);
          assert(branch != nullptr);
          if (geometry_t::intersects_box_box(min, max, branch->bmin(),
                                             branch->bmax())) {
            stk.push(branch);
          }
        }
      }
    }
    return ents;
  }

  void get_leaves(std::vector<branch_t *> &leaves) {
    for (auto &it : branch_map_) {
      if (it.second.is_leaf()) {
        leaves.push_back(&(it.second));
      }
    }
#ifdef DEBUG
    // Check if unique
    auto it = std::unique(leaves.begin(), leaves.end());
    assert(it == leaves.end());
#endif
  }

  void remove_non_local() {
    // remove the non local branches in the map
    auto it = branch_map_.begin();
    while (it != branch_map_.end()) {
      if (!it->second.is_local()) {
        auto parent_key = it->second.key();
        parent_key.pop();
        auto parent = branch_map_.find(parent_key);
        if (parent != branch_map_.end())
          parent->second.remove_bit(it->second.key().last_value());
        it = branch_map_.erase(it);
      } else {
        ++it;
      }
    }
  }

  /*!
    Construct a new entity. The entity's constructor should not be called
    directly.
   */
  template <class... Args> size_t make_entity(Args &&... args) {
    tree_entities_.emplace_back(std::forward<Args>(args)...);
    auto ent = &(tree_entities_.back());
    // Size -1 to start at 0
    size_t id = tree_entities_.size() - 1;
    ent->set_id_(id);
    return id;
  }

  /**
   * Insert directly a branch (certainly remote) in the tree
   */
  void insert_branch(const point_t &coordinates, const element_t &mass,
                     const point_t &bmin, const point_t &bmax, const key_t &key,
                     const int &owner, const size_t &sub_entities) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(rank != owner);
    // Check if this key already exists
    auto itr = branch_map_.find(key);
    // Case 1, branch does not exists localy
    if (itr == branch_map_.end()) {
      // Add the missing parents
      key_t pk = key;
      int last_bit = pk.last_value();
      pk.pop();
      while (branch_map_.find(pk) == branch_map_.end()) {
        branch_map_.emplace(pk, pk);
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
        itr->second.add_bit_child(last_bit);
        last_bit = pk.last_value();
        pk.pop();
      }
      // Set upper level not to leave
      branch_map_.find(pk)->second.set_leaf(false);
      branch_map_.find(pk)->second.add_bit_child(last_bit);

      branch_map_.emplace(key, key);
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

      size_t depth = key.depth();
      // Set the new depth of the tree
      max_depth_ = std::max(max_depth_,depth);
    } else {
      if (itr->second.owner() == rank)
        assert(itr->second.is_shared());
      else {
        // DO NOTHING, this branch have already been updated
        assert(!itr->second.is_local());
      }
    }
    // Add this branch if does not exists
  }

  /**
   * @brief Compute the keys of all the entities present in the structure
   */
  void compute_keys() {
    key_t::set_range(range_); 
#pragma omp parallel for
    for (size_t i = 0; i < entities_.size(); ++i) {
      entities_[i].set_key(key_t(entities_[i].coordinates()));
    }
  }

  /*!
    Return the tree's current max depth.
   */
  size_t max_depth() const { return max_depth_; }

  /*!
    Get an entity by entity id.
   */
  tree_entity_t *get(size_t id) {
    assert(id < tree_entities_.size());
    return &(tree_entities_[id]);
  }

  /**
   * @brief Get a branch by its id
   */
  branch_t *get(branch_id_t id) {
    auto itr = branch_map_.find(id);
    assert(itr != branch_map_.end());
    return &itr->second;
  }

  /*!
    Get the root branch (depth 0).
   */
  branch_t *root() { return &root_->second; }

  /**
   * @brief Generic information for the tree topology
   */
  friend std::ostream &operator<<(std::ostream &os, tree_topology &t) {
    os << "Tree: "
       << "#brchs: " << t.branch_map_.size()
       << " #ents: " << t.tree_entities_.size();
    os << " #root_subents: " << t.root()->sub_entities();
    os << " depth: " << t.max_depth_;
    return os;
  }

  /**
   * @brief      Export to a file the current tree in memory
   * This is useful for small number of particles to help representing the tree
   *
   * @param      tree   The tree to output
   * @param      range  The range of the particles, use to construct entity_key
   */
  void mpi_tree_traversal_graphviz(int num) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    clog_one(trace) << rank << " outputing tree file #" << num << std::endl;

    char fname[64];
    sprintf(fname, "output_graphviz_%02d_%02d.gv", rank, num);
    std::ofstream output;
    output.open(fname);
    output << "digraph G {" << std::endl << "forcelabels=true;" << std::endl;

    // Add the legend
    output << "branch [label=\"branch\" xlabel=\"sub_entities,owner\"]"
           << std::endl;

    std::stack<branch_t *> stk;
    // Get root
    auto rt = root();
    stk.push(rt);

    while (!stk.empty()) {
      branch_t *cur = stk.top();
      stk.pop();
      if (!cur->is_leaf()) {
        output << cur->key() << " [label=\"" << cur->key() << "\", xlabel=\""
               << cur->sub_entities() << " - " << cur->owner() << "\"];"
               << std::endl;
        switch (cur->locality()) {
        case 1:
          output << cur->key() << " [shape=circle,color=blue]" << std::endl;
          break;
        case 2:
          output << cur->key() << " [shape=circle,color=red]" << std::endl;
          break;
        case 3:
          output << cur->key() << " [shape=circle,color=green]" << std::endl;
          break;
        default:
          output << cur->key() << " [shape=circle,color=black]" << std::endl;
          break;
        }

        // Add the child to the stack and add for display
        for (size_t i = 0; i < (1 << dimension); ++i) {
          auto br = child(cur, i);
          if (br == nullptr)
            continue;
          stk.push(br);
          output << std::oct << cur->key() << "->" << br->key() << std::dec
                 << std::endl;
        }
      } else {
        output << cur->key() << " [label=\"" << cur->key() << "\", xlabel=\""
               << cur->sub_entities() << " - " << cur->owner() << "\"];"
               << std::endl;
        switch (cur->locality()) {
        case 1:
          output << cur->key() << " [shape=circle,color=blue]" << std::endl;
          break;
        case 2:
          output << cur->key() << " [shape=circle,color=red]" << std::endl;
          break;
        case 3:
          output << cur->key() << " [shape=circle,color=green]" << std::endl;
          break;
        default:
          output << cur->key() << " [shape=circle,color=black]" << std::endl;
          break;
        }
        for (auto ent : *cur) {
          auto e = get(ent);
          key_t key(range(), e->coordinates());
          key.truncate(max_depth() + 2);

          output << key << " [label=\"" << key << "\", xlabel=\"" << e->owner()
                 << " - " << e->global_id() << "\"];" << std::endl;

          output << cur->key() << "->" << key << std::endl;
          switch (e->locality()) {
          case 2:
            output << key << " [shape=box,color=green]" << std::endl;
            break;
          case 3:
            output << key << " [shape=box,color=black]" << std::endl;
            break;
          case 1:
            output << key << " [shape=box,color=red]" << std::endl;
            break;
          default:
            output << key << " [shape=circle,color=red]" << std::endl;
            break;
          }
          output << std::dec;
        }
      }
    }
    output << "}" << std::endl;
    output.close();
  }

  /**
   * @brief Try to insert an entity in the tree. This might need to refine
   * the branch.
   */
  void insert(const size_t &id) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // Find parent of the id
    auto ent = &(tree_entities_[id]);
    branch_id_t bid = ent->key();
    assert(bid.depth() > max_depth_);
    branch_t &b = find_parent(bid);
    // It is not a leaf, need to insert intermediate branch
    if (!b.is_leaf()) {
      // Create the branch
      size_t depth = b.key().depth() + 1;
      bid.truncate(depth);
      int bit = bid.last_value();
      b.add_bit_child(bit);
      branch_map_.emplace(bid, bid);
      branch_map_.find(bid)->second.set_leaf(true);
      branch_map_.find(bid)->second.insert(id);
    } else {
      // Conflict with a children
      if (b.size() == (1 << dimension)) {
        refine_(b);
        insert(id);
      } else {
        b.insert(id);
      }
    }
  }

private:
  /**
   * @brief Find the parent of an entity or branch based on the key
   * @details First truncate the key to the lowest possible in the tree, then
   * loop on the key to find an existing branch. At least it will find the root
   */
  branch_t &find_parent(branch_id_t bid) {
    branch_id_t pid = bid;
    pid.truncate(max_depth_);
    while (pid != root_->second.key()) {
      auto itr = branch_map_.find(pid);
      if (itr != branch_map_.end()) {
        return itr->second;
      }
      pid.pop();
    }
    return root_->second;
  }

  /**
   * @brief Refine the current branch b if there is a conflict of children
   */
  void refine_(branch_t &b) {
    branch_id_t pid = b.key();
    size_t depth = pid.depth() + 1;

    // For every children
    char bit_child = 0;
    for (auto ent : b) {
      key_t k = get(ent)->key();
      k.truncate(depth);
      bit_child |= 1 << k.last_value();
      branch_map_.emplace(k, k);
    }
    max_depth_ = std::max(max_depth_, depth);

    for (auto ent : b) {
      insert(ent);
    }

    b.set_leaf(false);
    b.clear();
    b.set_bit_child(bit_child);
  }

  // using branch_map_t = hashtable<key_int_t,branch_t>;
  using branch_map_t =
      std::unordered_map<branch_id_t, branch_t, branch_id_hasher__<key_t>>;

  branch_map_t branch_map_;
  size_t max_depth_;
  typename std::unordered_map<branch_id_t, branch_t,
                              branch_id_hasher__<key_t>>::iterator root_;
  // typename branch_map_t::iterator root_;
  range_t range_;
  std::vector<tree_entity_t> tree_entities_;
  std::vector<entity_t> entities_;
  std::vector<entity_t> entities_w_;

  const size_t max_traversal = 5;
  std::vector<std::vector<entity_t>> ghosts_entities_;
  size_t current_ghosts = 0;

  std::vector<entity_t> shared_entities_;

  const int ncritical = 32;
};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
