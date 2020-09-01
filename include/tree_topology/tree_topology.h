/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2019 Triad National Security, LLC
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

#ifndef flecsi_topology_tree_topology_h
#define flecsi_topology_tree_topology_h

/*!
  \file tree_topology.h
  \authors jloiseau@lanl.gov
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

#include "log.h"

#include "space_vector.h"

//#include "hashtable.h"
#include "tree_geometry.h"
#include "tree_types.h"

#ifdef ENABLE_DEBUG_TREE
#define _DEBUG_TREE_
#warning "Tree in debug mode with assert"
#endif

namespace flecsi {
namespace topology {

/*!
  The tree topology is parameterized on a policy P which defines its branch and
  entity types.
 */
template<class P>
class tree_topology : public P, public data::data_client_t
{

public:
  using Policy = P;

  static const size_t dimension = Policy::dimension;
  using element_t = typename Policy::element_t;
  using point_t = space_vector_u<element_t, dimension>;
  using range_t = std::array<point_t, 2>;
  using key_t = typename Policy::key_t;
  using entity_t = typename Policy::entity_t;
  using geometry_t = tree_geometry<element_t, dimension>;
  using cofm_t = typename Policy::cofm_t;
  using hcell_t = hcell<dimension, key_t, cofm_t, entity_t>;
  using key_int_t = typename Policy::key_int_t;

private:
  /**
   * @brief Entity type for MPI communication.
   * It requires the entity, its key in the tree (not full key) and
   * the rank that owns this entity for later communications.
   */
  struct share_entity_t {
    share_entity_t() {}
    share_entity_t(const int & o, const key_t & k, const entity_t & e)
      : owner(o), key(k), entity(e){};
    int owner;
    key_t key;
    entity_t entity;
  };
  /**
   * @brief Node type for MPI communication.
   * It requires the node (cofm), its keey in the tree and the rank
   * that owns this node.
   */
  struct share_node_t {
    share_node_t() {}
    share_node_t(const int & o, const key_t & k, const cofm_t & n, int nc)
      : owner(o), key(k), node(n), nchildren(nc) {};
    int owner;
    key_t key;
    cofm_t node;
    int nchildren;
  };

  /**
   * @brief Types for MPI communications
   * REQUEST: send a key request to another rank
   * REPLY_NODE: reply to another rank request with nodes
   * REPLY_ENTITY: reply to another rank request with entities
   * DONE_COMMS: Local rank done, send notification to other ranks
   */
  enum COMMS : int {
    REQUEST = 10,
    REQUEST_SUBTREE = 11,
    REPLY_NODE = 12,
    REPLY_ENTITY = 13,
    DONE_COMMS = 14
  };

public:
  tree_topology() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    comms_done_.resize(size);
  }
  ~tree_topology() {}

  /**
   * Clean the tree topology but not the local bodies
   * Remove shared entites and center of masses
   */
  void clean() {
    cofm_.clear();
    htable_.clear();
    shared_entities_.clear();
    shared_nodes_.clear();
  }

  /**
   * @brief Reset the ghosts, clean the tree and reconstruct it.
   * Do not share the particles again, use the current version of the keys
   */
  template<typename CCOFM>
  void reset_ghosts(CCOFM && f_c, bool do_share_edge = true) {
    clean();
    build_tree(f_c);
  }

  /**
   * \brief Change the range of the tree topology
   */
  void set_range(const range_t & range) {
    range_ = range;
  }

  /**
   * @brief Get the range
   */
  const std::array<point_t, 2> & range() {
    return range_;
  }

  /**
   * @ brief Return a reference to the vector of the entities
   */
  std::vector<entity_t> & entities() {
    return entities_;
  }

  /**
   * @brief Return an entity by its id
   */
  template<typename E>
  entity_t & entity(E e) {
    return entities_[static_cast<int>(e)];
  }

  /**
   * @brief Generic traversal function
   */
  template<typename FUNC, typename... ARGS>
  void traversal(hcell_t * cell, FUNC && func, ARGS &&... args) {
    std::stack<hcell_t *> stk;
    stk.push(cell);
    while(!stk.empty()) {
      hcell_t * cur = stk.top();
      stk.pop();
      if(func(cur, std::forward<ARGS>(args)...)) {
        hcell_t * daughters[nchildren_] = {nullptr};
        int children = 0;
        daughters_(cur, daughters, children);
        for(int i = 0; i < children; ++i) {
          stk.push(daughters[i]);
        } // for
      } // if
    } // while
  }

  /**
` * @brief Apply a function EF to the sub_cells using asynchronous comms.
  */
  template<typename EF, typename... ARGS>
  void traversal_sph(EF && ef, ARGS &&... args) {
    log_one(trace) << "Traversal SPH" << std::endl;
    double start = omp_get_wtime();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Find all nodes of the tree with at most sub_entities_ elements
    std::vector<key_t> cells;
    traversal(
      root(),
      [&](hcell_t * cell, std::vector<key_t> & c, const int & sent) {
        if(cell->is_node() &&
           (cell->is_shared() || get_node(cell)->sub_entities() > sent)) {
          return true;
        }
        if(cell->is_node() || (cell->is_entity() && !cell->is_shared())) {
          c.push_back(cell->key());
        }
        return false;
      } // lambda
      ,
      cells, sub_entities_);

    // prepare comms arrays
    init_comms_(size);
    std::stack<key_t> stk_nonlocal;

    // Traversal data
    std::vector<std::vector<key_t>> request_keys;
    request_keys.resize(size);
    std::vector<hcell_t *> * queue = new std::vector<hcell_t *>();
    std::vector<hcell_t *> * new_queue = new std::vector<hcell_t *>();
    std::vector<std::vector<entity_t *>> neighbors;
    hcell_t* daughters[nchildren_];
    int children;

    int i = 0;
    double lost_time;
    bool alternate = true;
    while(i < cells.size() || !stk_nonlocal.empty()) {

      if(size > 1)
        check_comms_();

      key_t curkey = key_t(0);
#ifdef _DEBUG_TREE_
      lost_time = omp_get_wtime();
#endif
      if(i >= cells.size())
        alternate = false;
      if(alternate) {
        curkey = cells[i++];
        alternate = false;
      }
      else {
        if(!stk_nonlocal.empty()) {
          curkey = stk_nonlocal.top();
          stk_nonlocal.pop();
        }
        else {
          if(i < cells.size())
            curkey = cells[i++];
          else
            break;
        }
        alternate = true;
      } // if
#ifdef _DEBUG_TREE_
      assert(curkey != key_t(0));
#endif
      bool non_local = false;
      bool rank_request = false;

      hcell_t * cur = &(htable_.find(curkey)->second);
      std::vector<entity_t *> cur_entities;

      cofm_t * cur_node = nullptr;

      if(cur->is_node()) {
        traversal(
          cur,
          [&](hcell_t * cell, std::vector<entity_t *> & ce) {
            if(cell->is_node()) {
              return true;
            }
            else {
              if(!cell->is_shared())
                ce.push_back(get_entity(cell));
            }
            return false;
          },
          cur_entities); // lambda
        cur_node = get_node(cur);
      }
      else {
        cur_entities.push_back(get_entity(cur));
      } // if

      neighbors.clear();
      neighbors.resize(cur_entities.size());
      queue->clear();
      queue->push_back(root());

      while(!queue->empty()) {
        new_queue->clear();
        // Eliminate geometrically
        for(int j = 0; j < queue->size(); ++j) {
          bool accepted = false;
          hcell_t * hcur = (*queue)[j];
          if(hcur->is_node()) {
            cofm_t * c = get_node(hcur);
            // Check if node concerned
            if(cur_node != nullptr) {
              if(!geometry_t::intersects_box_box(
                   c->bmin(), c->bmax(), cur_node->bmin(), cur_node->bmax())) {
                continue;
              }
            } // if
            // If yes, check for all entities before request
            for(int k = 0; k < cur_entities.size() && !accepted; ++k) {
              if(geometry_t::intersects_sphere_box(c->bmin(), c->bmax(),
                   cur_entities[k]->coordinates(), cur_entities[k]->radius())) {
                accepted = true;
                if(hcur->is_empty_node()) {
                  non_local = true;
                  if(!hcur->requested()) {
#ifdef _DEBUG_TREE_
                    assert(hcur->owner() != rank);
#endif
                    hcur->set_requested();
                    request_keys[hcur->owner()].push_back(hcur->key());
                    rank_request = true;
                  }
                }
                else {
                  children = 0;
                  daughters_(hcur, daughters, children);
                  for(int l = 0; l < children; ++l)
                    new_queue->push_back(daughters[l]);
                } // if
              } // if
            } // if
          }
          else {
#ifdef _DEBUG_TREE_
            assert(hcur->is_entity());
#endif
            entity_t * e = get_entity(hcur);
#ifdef _DEBUG_TREE_
            assert(e != nullptr);
#endif
            if(cur_node != nullptr) {
              element_t extent_ent =
                std::max(e->radius(), cur_node->lap()) + cur_node->radius();
              if(!geometry_t::within_distance2(
                   e->coordinates(), cur_node->coordinates(), extent_ent))
                continue;
            }
            for(int k = 0; k < cur_entities.size(); ++k) {
              element_t extent =
                std::max(cur_entities[k]->radius(), e->radius());
              if(geometry_t::within_distance2(
                   cur_entities[k]->coordinates(), e->coordinates(), extent)) {
                neighbors[k].push_back(e);
              } // if
            } // for
          } // if
        } // for
        if(non_local) {
          if(rank_request) {
            request_(request_keys);
            for(int k = 0; k < size; ++k) {
              request_keys[k].clear();
            } // for
          } // if
#ifdef _DEBUG_TREE_
          lost_timer_ += omp_get_wtime() - lost_time;
#endif
          stk_nonlocal.push(curkey);
          break;
        } // if

        auto tmp = queue;
        queue = new_queue;
        new_queue = tmp;

      } // while
      if(!non_local) {
        for(int j = 0; j < cur_entities.size(); ++j) {
#ifdef _DEBUG_TREE_
          assert(neighbors[j].size() != 0);
#endif
          ef(*cur_entities[j], neighbors[j], std::forward<ARGS>(args)...);
        } // for
      } // if
    } // while
    if(size > 1) {
      comms_all_done_ = false;
      std::vector<MPI_Request> done_requests(size);
      std::vector<MPI_Status>  done_status(size);
      for(int i = 0; i < size; ++i) {
        MPI_Issend(nullptr, 0, MPI_INT, i, DONE_COMMS, MPI_COMM_WORLD,
            &done_requests[i]);
      } // for
      while(!comms_all_done_) {
        wait_comms_();
      } // while
      MPI_Waitall(size, &done_requests[0], &done_status[0]);
    } // if

    clean_comms_();

    queue->clear();
    new_queue->clear();
    delete queue;
    delete new_queue;

    MPI_Barrier(MPI_COMM_WORLD);
    double tree_timer = omp_get_wtime() - start;
    log_one(trace) << std::fixed << std::setprecision(3)
                   << "Traversal SPH.done: " << tree_timer << "s"
#ifdef _DEBUG_TREE_
                   << " comms_: " << comms_timer_ << "s ("
                   << comms_timer_ * 100 / tree_timer << "%) "
                   << "lost_: " << lost_timer_ << "s ("
                   << lost_timer_ * 100 / tree_timer << "%)"
#endif
                   << std::endl;
  } // traversal_sph

  /**
   * @brief Fast Multipole Method Traversal.
   * Perform a tree traversal and update the missing neighbors.
   */
  template<typename C2C, typename P2C, typename P2P, typename C2P>
  void traversal_fmm(const double MAC,
    C2C && t_c2c,
    P2C && t_p2c,
    P2P && f_p2p,
    C2P && f_c2p) {
    log_one(trace) << "Traversal FMM (" << MAC << ")" << std::endl;
    double start = omp_get_wtime();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    init_comms_(size);

    // Find pairs of interacting cells
    using interaction_t = std::pair<key_t, key_t>;
    std::vector<interaction_t> * queue = new std::vector<interaction_t>();
    std::vector<interaction_t> * new_queue = new std::vector<interaction_t>();
    std::vector<interaction_t> p2p;
    std::vector<entity_t *> subs;
    std::vector<entity_t *> neighbors;
    hcell_t * daughters[nchildren_];
    int children;
    double lost_time;

    std::vector<std::vector<key_t>> request_keys;
    request_keys.resize(size);

    queue->emplace_back(key_t::root(), key_t::root());
    while(not queue->empty()) {

      if(size > 1)
        check_comms_();

      bool rank_request = false;

      new_queue->clear();
      for(int i = 0; i < queue->size(); ++i) {

#ifdef _DEBUG_TREE_
        lost_time = omp_get_wtime();
#endif

        key_t khc1 = (*queue)[i].first;
        key_t khc2 = (*queue)[i].second;
        hcell_t * hc1 = &(htable_.find(khc1)->second);
        hcell_t * hc2 = &(htable_.find(khc2)->second);

        assert(hc1->iam_owner());

        if(!hc2->is_empty_node()) {
          if(hc1->is_entity() && hc2->is_entity()) {
            // both are entities: append interaction to the p2p list
            p2p.push_back((*queue)[i]);
          }
          else { // at least one is a node

            if(hc1->key() == hc2->key()) { // same node
              // check for the number of subentities

              if(get_node(hc1)->sub_entities() < fmm_sub_entities_) {
                p2p.push_back((*queue)[i]);
              }
              else {
                // split it for self-interaction
                daughters_(hc1, daughters, children);
                for(int k1 = 0; k1 < children; ++k1) {
                  if(daughters[k1]->iam_owner())
                    new_queue->emplace_back(
                      daughters[k1]->key(), daughters[k1]->key());
                  for(int k2 = k1 + 1; k2 < children; ++k2) {
                    if(daughters[k1]->iam_owner())
                      new_queue->emplace_back(
                        daughters[k1]->key(), daughters[k2]->key());
                    if(daughters[k2]->iam_owner())
                      new_queue->emplace_back(
                        daughters[k2]->key(), daughters[k1]->key());
                  }
                } // for k1
              }
            }
            else { // different nodes
              point_t coords1 = {};
              element_t radius1 = 0;
              int subent1 = 1;
              if(hc1->is_node()) {
                cofm_t * n = get_node(hc1);
                coords1 = n->coordinates();
                radius1 = n->radius();
                subent1 = n->sub_entities();
              }
              else {
                entity_t * e = get_entity(hc1);
                coords1 = e->coordinates();
              }

              point_t coords2 = {};
              element_t radius2 = 0;
              int subent2 = 1;
              if(hc2->is_node()) {
                cofm_t * n = get_node(hc2);
                coords2 = n->coordinates();
                radius2 = n->radius();
                subent2 = n->sub_entities();
              }
              else {
                entity_t * e = get_entity(hc2);
                coords2 = e->coordinates();
              }

              if(geometry_t::mac(coords1, radius1, coords2, radius2, MAC)) {
                assert(hc1->is_node() or hc2->is_node());
                if(hc1->is_node()) {
                  cofm_t * n1 = get_node(hc1);
                  if(hc2->is_node()) {
                    t_c2c(n1, get_node(hc2));
                  }
                  else {
                    entity_t * e = get_entity(hc2);
                    t_p2c(n1, e);
                  }
                  // save this node for later c2c interactions
                  n1->set_affected(true);
                }
                else { // hc1 is an entity
                  neighbors.clear();
                  subs.clear();
                  subs.push_back(get_entity(hc1));
                  f_p2p(subs, get_node(hc2), neighbors);
                }
              }
              else { // nodes do not satisfy MAC
                if(subent1 + subent2 < fmm_sub_entities_) {
                  // if not enough subentities, give up with splitting
                  p2p.push_back((*queue)[i]);
                  std::vector<std::vector<key_t>> request_keys_subtree(size);
                  bool rqst_subtree = false;
                  if(hc2->is_shared()) {
                    traversal(
                      hc2,
                      [&](
                        hcell_t * cell, std::vector<std::vector<key_t>> & nk) {
                        if((cell->is_node() && !cell->is_shared()) ||
                           cell->is_entity()) {
                          return false;
                        }
                        // if(cell->is_node() && cell->is_shared()){
                        //  return true;
                        //}
                        if(cell->is_empty_node() && !cell->requested()) {
                          rqst_subtree = true;
                          assert(cell->owner() != rank);
                          cell->set_requested();
                          nk[cell->owner()].push_back(cell->key());
                          return false;
                        }
                        return true;
                      } // lambda
                      ,
                      request_keys_subtree);
                  }
                  // Send request
                  if(rqst_subtree)
                    request_(request_keys_subtree, REQUEST_SUBTREE);
                  // Retrieve the non local particles of this sub-tree
                }
                else {
                  if(radius1 > radius2) { // split the bigger node
                    // node that if one of the cells is an entity, then its
                    // radius will be zero; the other one must be the node with
                    // nonzero radius
                    daughters_(hc1, daughters, children);
                    for(int k = 0; k < children; ++k) {
                      if(daughters[k]->iam_owner()) {
                        new_queue->emplace_back(
                          daughters[k]->key(), hc2->key());
                      }
                    }
                  }
                  else {
                    daughters_(hc2, daughters, children);
                    for(int k = 0; k < children; ++k) {
                      new_queue->emplace_back(hc1->key(), daughters[k]->key());
                    }
                  }
                } // if enough subentities for splitting
              } // if not MAC
            } // if different nodes
          } // if at least one is a node
        }
        else {
          // Check if node is empty and retrieve if needed
          if(!hc2->requested()) {
#ifdef _DEBUG_TREE_
            assert(hc2->owner() != rank);
#endif
            hc2->set_requested();
            request_keys[hc2->owner()].push_back(hc2->key());
            rank_request = true;
          }
          new_queue->emplace_back(hc1->key(), hc2->key());
#ifdef _DEBUG_TREE_
          lost_timer_ += omp_get_wtime() - lost_time;
#endif
        } // if
      } // loop over the queue
      if(rank_request) {
        request_(request_keys);
        for(int k = 0; k < request_keys.size(); ++k) {
          request_keys[k].clear();
        }
      } // if non_local
      auto tmp = queue;
      queue = new_queue;
      new_queue = tmp;
    } // while queue

    if(size > 1) {
      comms_all_done_ = false;
      std::vector<MPI_Request> done_requests(size);
      std::vector<MPI_Status>  done_status(size);
      for(int i = 0; i < size; ++i) {
        MPI_Issend(nullptr, 0, MPI_INT, i, DONE_COMMS, MPI_COMM_WORLD,
            &done_requests[i]);
      } // for
      // Handle communications
      while(!comms_all_done_) {
        wait_comms_();
      } // while
      MPI_Waitall(size, &done_requests[0], &done_status[0]);
    }


    // node-node interaction
    std::vector<hcell_t *> affected_nodes;
    traversal(
      root(),
      [&](hcell_t * cell, std::vector<hcell_t *> & hc) {
        if(!cell->iam_owner()) {
          return false; // do not expand others' nodes
        }
        if(cell->is_node() && get_node(cell)->affected()) {
          hc.push_back(cell);
        }
        return true;
      } // lambda
      ,
      affected_nodes);

    neighbors.clear();
    for(int i = 0; i < affected_nodes.size(); ++i) {
      hcell_t * hc = affected_nodes[i];
      subs.clear();

      // Find all sub entities
      traversal(
        hc,
        [&](hcell_t * cell, std::vector<entity_t *> & e) {
          if(cell->is_node()) {
            return true;
          }
          if(cell->is_entity() && !cell->is_shared()) {
            e.push_back(get_entity(cell));
          }
          return false;
        } // lambda
        ,
        subs);

      f_c2p(get_node(hc), subs);
    }

    for(int i = 0; i < p2p.size(); ++i) {
      hcell_t * hc1 = &(htable_.find(p2p[i].first)->second);
      hcell_t * hc2 = &(htable_.find(p2p[i].second)->second);

      // subentities of hc1
      std::vector<entity_t *> subs;
      if(hc1->is_node()) {
        traversal(
          hc1,
          [&](hcell_t * cell, std::vector<entity_t *> & e) {
            if(cell->is_node()) {
              return true;
            }
            if(cell->is_entity() && !cell->is_shared()) {
              e.push_back(get_entity(cell));
            }
            return false;
          } // lambda
          ,
          subs);
      }
      else {
        subs.push_back(get_entity(hc1));
      }

      // use 'neighbors' vector to store subentities of hc2
      neighbors.clear();
      if(hc2->is_node()) {
        traversal(
          hc2,
          [&](hcell_t * cell, std::vector<entity_t *> & e) {
            if(cell->is_node()) {
              return true;
            }
            if(cell->is_entity() && !cell->is_shared()) {
              e.push_back(get_entity(cell));
            }
            return false;
          } // lambda
          ,
          neighbors);
      }
      else {
        neighbors.push_back(get_entity(hc2));
      }

      if(hc1->is_node()) {
        f_c2p(get_node(hc1), subs);
      }
      else {
        subs.clear();
        subs.push_back(get_entity(hc1));
        f_p2p(subs, nullptr, neighbors);
      }

    } // for p2p interactions

    clean_comms_();

    queue->clear();
    new_queue->clear();
    delete queue;
    delete new_queue;

    MPI_Barrier(MPI_COMM_WORLD);
    double tree_timer = omp_get_wtime() - start;
    log_one(trace) << std::fixed << std::setprecision(3)
                   << "Traversal FMM.done: " << tree_timer << "s"
#ifdef _DEBUG_TREE_
                   << " comms_: " << comms_timer_ << "s ("
                   << comms_timer_ * 100 / tree_timer << "%) "
                   << "lost_: " << lost_timer_ << "s ("
                   << lost_timer_ * 100 / tree_timer << "%)"
#endif
                   << std::endl;
  }

  /**
   * @brief return a vector of entities in the specified spheroid
   */
  template<typename EF>
  std::vector<entity_t *>
  find_in_radius(const point_t & center, element_t radius, EF && ef) {
    std::vector<entity_t *> result;
    traversal(
      root(),
      [&](hcell_t * cur, std::vector<entity_t *> & result) {
        if(cur->is_node()) {
          cofm_t * c = get_node(cur);
          element_t extent = std::max(c->lap(), radius) + c->radius();
          if(geometry_t::within_distance2(c->coordinates(), center, extent))
            return true;
        }
        else {
          entity_t * e = get_entity(cur);
          element_t extent = std::max(radius, e->radius());
          if(geometry_t::within_distance2(center, e->coordinates(), extent))
            result.push_back(e);
        }
        return false;
      },
      result);
    return result;
  }

  /**
   * @brief Compute the keys of all the entities present in the structure
   */
  void compute_keys() {
    for(size_t i = 0; i < entities_.size(); ++i) {
      entities_[i].set_key(key_t(range_, entities_[i].coordinates()));
    } // for
  }

  /*!
    @brief eturn the tree's current max depth.
   */
  size_t max_depth() const {
    return max_depth_;
  }

  /*!
    @brief Get the root branch (depth 0).
   */
  hcell_t * root() {
    return &root_->second;
  }

  cofm_t * root_node() {
    return get_node(root());
  }

  /**
   * @brief Generic information for the tree topology
   */
  friend std::ostream & operator<<(std::ostream & os, tree_topology & t) {
    auto r = t.htable_.find(key_t::root());
    cofm_t * root_ptr = r->second.is_shared()
                          ? &t.shared_nodes_[r->second.node_idx()]
                          : &(t.cofm_[r->second.node_idx()]);
    os << "Tree: "
       << "#node: " << t.htable_.size() - t.entities_.size();
    os << " depth: " << t.max_depth_;
    os << " #root_subents: " << root_ptr->sub_entities();
    // os << " center: "<<root_ptr->coordinates();
    os << " mass: " << root_ptr->mass();
    os << " radius: " << root_ptr->radius();
    return os;
  }

  /**
   * @brief Loop over the bodies to insert them in the tree and construct the
   * nodes.
   * 1. Exchange the boundaries with the neighbors
   * 2. insert the entities and create the branches
   * 2.a. If a branch is between lo-hi key, the cofm can be computed
   * 3. The tree is ready to share entities/nodes with the neighbors
   **/
  template<typename CCOFM>
  void build_tree(CCOFM && f_cc) {
    log_one(trace) << "Building tree" << std::endl;
    double start = omp_get_wtime();
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Exchange high and low bound */
    key_t lokey = entities_[0].key();
    key_t hikey = entities_[entities_.size() - 1].key();
    exchange_boundaries_(hikey, lokey, hibound_, lobound_);
    max_depth_ = 0;
    // Add the root
    htable_.emplace(key_t::root(), key_t::root());
    root_ = htable_.find(key_t::root());

    size_t current_depth = key_t::max_depth();
    // Entity keys, last and current
    key_t lastekey = key_t(0);
    if(rank != 0)
      lastekey = lobound_;
    key_t ekey;
    // Node keys, last and Current
    key_t lastnkey = key_t::root();
    key_t nkey, loboundnode, hiboundnode;
    // Current parent and value
    hcell_t * parent = nullptr;
    int oldidx = -1;

    bool iam0 = rank == 0;
    bool iamlast = rank == size - 1;

#ifdef _DEBUG_TREE_
    assert(lobound_ <= lokey);
    assert(hibound_ >= hikey);
#endif

    // The extra turn in the loop is to finish the missing
    // parent of the last entity
    for(size_t i = 0; i <= entities_.size(); ++i) {
      if(i < entities_.size()) {
        ekey = entities_[i].key();
        // Compute the current node key
      }
      else {
        ekey = hibound_;
      }
      nkey = ekey;
      nkey.pop(current_depth);
      bool loopagain = false;
      // Loop while there is a difference in the current keys
      while(nkey != lastnkey || (iamlast && i == entities_.size())) {
        loboundnode = lobound_;
        loboundnode.pop(current_depth);
        hiboundnode = hibound_;
        hiboundnode.pop(current_depth);
        if(loopagain && (iam0 || lastnkey > loboundnode) &&
           (iamlast || lastnkey < hiboundnode)) {
          // This node is done, we can compute CoFM
          finish_(lastnkey, f_cc);
        }
        if(iamlast && lastnkey == key_t::root())
          break;
        loopagain = true;
        current_depth++;
        nkey = ekey;
        nkey.pop(current_depth);
        lastnkey = lastekey;
        lastnkey.pop(current_depth);
      } // while

      if(iamlast && i == entities_.size())
        break;

      parent = &(htable_.find(lastnkey)->second);
      oldidx = parent->entity_idx();
      // Insert the eventual missing parents in the tree
      // Find the current parent of the two entities
      while(1) {
        current_depth--;
        lastnkey = lastekey;
        lastnkey.pop(current_depth);
        nkey = ekey;
        nkey.pop(current_depth);
        if(nkey != lastnkey)
          break;
        // Add a children
        int bit = nkey.last_value();
        parent->add_child(bit);
        parent->set_entity_idx(-1);
        htable_.emplace(nkey, nkey);
        parent = &(htable_.find(nkey)->second);
      } // while

      // Recover deleted entity
      if(oldidx != -1) {
        int bit = lastnkey.last_value();
        parent->add_child(bit);
        parent->set_entity_idx(-1);
        htable_.emplace(lastnkey, hcell_t(lastnkey, i - 1));
      } // if

      if(i < entities_.size()) {
        // Insert the new entity
        int bit = nkey.last_value();
        parent->add_child(bit);
        htable_.emplace(nkey, hcell_t(nkey, i));
      } // if

      // Prepare next loop
      lastekey = ekey;
      lastnkey = nkey;
      max_depth_ = std::max(max_depth_, current_depth);
    } // for
    share_nodes_(f_cc);
    MPI_Barrier(MPI_COMM_WORLD);
    log_one(trace) << "Building tree.done: " << omp_get_wtime() - start << "s"
                   << std::endl;
  }

  /**
   * @brief Return an entity linked to a cell
   * This takes care of the local/shared entity
   */
  entity_t * get_entity(const hcell_t * hc) {
#ifdef _DEBUG_TREE_
    assert(hc->is_entity());
#endif
    int idx = hc->entity_idx();
#ifdef _DEBUG_TREE_
    assert(
      hc->is_shared() ? idx < shared_entities_.size() : idx < entities_.size());
#endif
    return hc->is_shared() ? &shared_entities_[idx] : &entities_[idx];
  }

  /**
   * @brief Return a node linked to a cell
   * This takes care of the local/shared node
   */
  cofm_t * get_node(const hcell_t * hc) {
#ifdef _DEBUG_TREE_
    assert(hc->is_node());
#endif
    int idx = hc->node_idx();
#ifdef _DEBUG_TREE_
    assert(hc->is_shared() ? idx < shared_nodes_.size() : idx < cofm_.size());
#endif
    if(hc->is_shared())
      return &shared_nodes_[idx];
    else
      return &cofm_[idx];
    assert(false);
    return nullptr;
  }

private:
  /**
   * @brief      Export to a file the current tree in memory
   * This is useful for small number of particles to see the tree
   * representation
   */
  void graphviz_draw_(int num) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    log_one(trace) << rank << " outputing tree file #" << num << std::endl;

    char fname[64];
    sprintf(fname, "output_graphviz_%02d_%02d.gv", rank, num);
    std::ofstream output;
    output.open(fname);
    output << "digraph G {" << std::endl << "forcelabels=true;" << std::endl;

    // Add the legend
    // output << "branch [label=\"branch\" xlabel=\"sub_entities,owner\"]"
    //       << std::endl;

    std::stack<hcell_t *> stk;
    // Get root
    stk.push(root());

    while(!stk.empty()) {
      hcell_t * cur = stk.top();
      stk.pop();
      if(cur->is_unset()) {
        output << std::oct << cur->key() << std::dec << " [label=\"" << std::oct
               << cur->key() << std::dec << "\", xlabel=\"\"];" << std::endl;
        output << std::oct << cur->key() << std::dec
               << " [shape=circle,color=black]" << std::endl;

        // Add the child to the stack and add for display
        for(size_t i = 0; i < nchildren_; ++i) {
          if(cur->get_child(i)) {
            key_t ckey = cur->key();
            ckey.push(i);
            auto it = htable_.find(ckey);
            if(it != htable_.end()) {
              stk.push(&it->second);
            }
            else {
              continue;
            }
            output << std::oct << cur->key() << "->" << it->second.key()
                   << std::dec << std::endl;
          } // if
        }
      }
      else if(cur->is_node()) {
        int sub_ent = 0;
        int idx = cur->node_idx();
        cofm_t * c = cur->is_shared() ? &shared_nodes_[idx] : &cofm_[idx];
        output << std::oct << cur->key() << std::dec << " [label=\"" << std::oct
               << cur->key() << std::dec << "\", xlabel=\"" << cur->nchildren()
               << "," << c->sub_entities() << "," << cur->owner() << "\"];"
               << std::endl;
        if(cur->is_shared()) {
          output << std::oct << cur->key() << std::dec
                 << " [shape=circle,color=green]" << std::endl;
        }
        else {
          output << std::oct << cur->key() << std::dec
                 << " [shape=circle,color=blue]" << std::endl;
        }

        // Add the child to the stack and add for display
        for(size_t i = 0; i < nchildren_; ++i) {
          if(cur->get_child(i)) {
            key_t ckey = cur->key();
            ckey.push(i);
            auto it = htable_.find(ckey);
            if(it != htable_.end()) {
              stk.push(&it->second);
            }
            else {
              continue;
            }
            output << std::oct << cur->key() << "->" << it->second.key()
                   << std::dec << std::endl;
          } // if
        }
      }
      else {
        output << std::oct << cur->key() << std::dec << " [label=\"" << std::oct
               << cur->key() << std::dec << "\", xlabel=\"" << cur->owner()
               << "\"];" << std::endl;
        if(cur->is_shared()) {
          output << std::oct << cur->key() << std::dec
                 << " [shape=circle,color=grey]" << std::endl;
        }
        else {
          output << std::oct << cur->key() << std::dec
                 << " [shape=circle,color=red]" << std::endl;
        }
      } // if
    } // while
    output << "}" << std::endl;
    output.close();
  }

  /**
   * @brief Check for communciation: requests or replies from other
   * ranks.
   */
  void check_comms_() {
#ifdef _DEBUG_TREE_
    double start = omp_get_wtime();
#endif
    int flag = 1, size, rank;
    MPI_Status status;
    // static int tree_num = 1 ;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // bool updated_tree = false;
    // Handle all current requests
    while(flag == 1) {
      // Change to MPI_Probe when replying only
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
      if(flag) {
        int source = status.MPI_SOURCE;
        int tag = status.MPI_TAG;
        int nrecv = 0;
#ifdef _DEBUG_TREE_
        if(tag != DONE_COMMS)
          assert(source != rank);
#endif
        MPI_Get_count(&status, MPI_BYTE, &nrecv);
        switch(tag) {
          case REQUEST_SUBTREE:
            recv_requests_subtree_(source, nrecv);
            break;
          case REQUEST:
            recv_requests_(source, nrecv);
            break;
          case REPLY_NODE:
            // updated_tree = true;
            recv_node_replies_(source, nrecv);
            break;
          case REPLY_ENTITY:
            // updated_tree = true;
            recv_entity_replies_(source, nrecv);
            break;
          case DONE_COMMS:
            MPI_Recv(nullptr, 0, MPI_INT, source, DONE_COMMS, MPI_COMM_WORLD,
              MPI_STATUS_IGNORE);
            comms_done_[source] = true;
            comms_all_done_ = true;
            for(int i = 0; i < size; ++i) {
              if(!comms_done_[i]) {
                comms_all_done_ = false;
                break;
              } // if
            } // for
            break;
          default:
            std::cerr << "Unknown message type: " << tag
                      << " source: " << source << std::endl;
            MPI_Finalize();
            exit(1);
        } // switch
      } // if
    } // while
#ifdef _DEBUG_TREE_
    comms_timer_ += omp_get_wtime() - start;
#endif
    // if(updated_tree){
    //  graphviz_draw(tree_num++);
    //}
  }

  void wait_comms_() {
#ifdef _DEBUG_TREE_
    double start = omp_get_wtime();
#endif
    int size, rank;
    bool end = false;
    MPI_Status status;
    // static int tree_num = 1 ;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // bool updated_tree = false;
    // Handle all current requests
    while(!comms_all_done_) {
      // Change to MPI_Probe when replying only
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      int source = status.MPI_SOURCE;
      int tag = status.MPI_TAG;
      int nrecv = 0;
#ifdef _DEBUG_TREE_
      if(tag != DONE_COMMS)
        assert(source != rank);
#endif
      MPI_Get_count(&status, MPI_BYTE, &nrecv);
      switch(tag) {
        case REQUEST_SUBTREE:
          recv_requests_subtree_(source, nrecv);
          break;
        case REQUEST:
          recv_requests_(source, nrecv);
          break;
        case REPLY_NODE:
          // updated_tree = true;
          recv_node_replies_(source, nrecv);
          break;
        case REPLY_ENTITY:
          // updated_tree = true;
          recv_entity_replies_(source, nrecv);
          break;
        case DONE_COMMS:
          MPI_Recv(nullptr, 0, MPI_INT, source, DONE_COMMS, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);
          comms_done_[source] = true;
          comms_all_done_ = true;
          for(int i = 0; i < size; ++i) {
            if(!comms_done_[i]) {
              comms_all_done_ = false;
              break;
            } // if
          } // for
          break;
        default:
          std::cerr << "Unknown message type: " << tag << " source: " << source
                    << std::endl;
          MPI_Finalize();
          exit(1);
      } // switch
    } // while
#ifdef _DEBUG_TREE_
    comms_timer_ += omp_get_wtime() - start;
#endif
    // if(updated_tree){
    //  graphviz_draw(tree_num++);
    //}
  }

  /**
   * @brief Request a set of keys from another rank.
   */
  void request_(const std::vector<std::vector<key_t>> & keys,
    int rtype = REQUEST) {
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for(int i = 0; i < size; ++i) {
      int ksize = keys[i].size();
      if(ksize > 0) {
        if(mpi_requests_[current_requests_].size() + 1 >=
           requests_keys_max_ - 1) {
          current_requests_++;
          mpi_requests_.resize(current_requests_ + 1);
          mpi_requests_[current_requests_].reserve(requests_keys_max_);
        } // if
        requests_keys_.push_back(keys[i]);
        int cksize = requests_keys_.back().size();
        mpi_requests_[current_requests_].push_back(MPI_Request{});
        MPI_Issend(&requests_keys_.back()[0], ksize * sizeof(key_t), MPI_BYTE, i,
          rtype, MPI_COMM_WORLD, &mpi_requests_[current_requests_].back());
      } // if
    } // for

  }

  /**
   * @brief Check if another rank requested a group of keys.
   * In this version the reply will be split between the nodes
   * and the entities present in the requested keys.
   **/
  void recv_requests_(const int & partner, const int & nrecv) {
    bool found = false;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nkeys = nrecv / sizeof(key_t);
    std::vector<key_t> keys(nkeys);
    MPI_Recv(&keys[0], nrecv, MPI_BYTE, partner, REQUEST, MPI_COMM_WORLD,
      MPI_STATUS_IGNORE);
    std::vector<share_node_t> tmp_nodes_replies;
    std::vector<share_entity_t> tmp_entities_replies;
    for(int i = 0; i < nkeys; ++i) {
      hcell_t * cur = &(htable_.find(keys[i])->second);
#ifdef _DEBUG_TREE_
      assert(cur->is_node());
#endif
      tmp_nodes_replies.emplace_back(cur->owner(),cur->key(),*get_node(cur),
          cur->nchildren());
      for(int j = 0; j < nchildren_; ++j) {
        if(cur->get_child(j)) {
          key_t ckey = cur->key();
          ckey.push(j);
          auto child = htable_.find(ckey);
#ifdef _DEBUG_TREE_
          assert(child != htable_.end());
#endif
          if(child->second.is_node()) {
            tmp_nodes_replies.emplace_back(child->second.owner(),
              child->second.key(), *get_node(&child->second),
              child->second.nchildren());
          }
          else if(child->second.is_entity()) {
            tmp_entities_replies.emplace_back(child->second.owner(),
              child->second.key(), *get_entity(&child->second));
          }
#ifdef _DEBUG_TREE_
          else {
            assert(false);
          } // if
#endif
        } // if
      } // for
    } // for
    if(tmp_nodes_replies.size() != 0) {
      mpi_replies_[current_replies_].push_back(MPI_Request{});
      nodes_replies_.push_back(tmp_nodes_replies);
      MPI_Issend(&nodes_replies_[nodes_replies_.size() - 1][0],
        sizeof(share_node_t) * tmp_nodes_replies.size(), MPI_BYTE, partner,
        REPLY_NODE, MPI_COMM_WORLD, &mpi_replies_[current_replies_].back());
      found = true;
      if(mpi_replies_[current_replies_].size() >= requests_keys_max_ - 1) {
        current_replies_++;
        mpi_replies_.resize(current_replies_ + 1);
        mpi_replies_[current_replies_].reserve(requests_keys_max_);
      } // if
    } // if
    if(tmp_entities_replies.size() != 0) {
      mpi_replies_[current_replies_].push_back(MPI_Request{});
      entities_replies_.push_back(tmp_entities_replies);
      MPI_Issend(&entities_replies_[entities_replies_.size() - 1][0],
        sizeof(share_entity_t) * tmp_entities_replies.size(), MPI_BYTE, partner,
        REPLY_ENTITY, MPI_COMM_WORLD, &mpi_replies_[current_replies_].back());
      found = true;
      if(mpi_replies_[current_replies_].size() >= requests_keys_max_ - 1) {
        current_replies_++;
        mpi_replies_.resize(current_replies_ + 1);
        mpi_replies_[current_replies_].reserve(requests_keys_max_);
      } // if
    } // if
#ifdef _DEBUG_TREE_
    assert(found);
#endif
  }

  /**
   * @brief Check if another rank requested a group of keys.
   * In this version the reply will be split between the nodes
   * and the entities present in the requested keys.
   **/
  void recv_requests_subtree_(const int & partner, const int & nrecv) {
    bool found = false;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nkeys = nrecv / sizeof(key_t);
    std::vector<key_t> keys(nkeys);
    MPI_Recv(&keys[0], nrecv, MPI_BYTE, partner, REQUEST_SUBTREE,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::vector<share_node_t> tmp_nodes_replies;
    std::vector<share_entity_t> tmp_entities_replies;
    for(int i = 0; i < nkeys; ++i) {
      hcell_t * cur = &(htable_.find(keys[i])->second);
#ifdef _DEBUG_TREE_
      assert(cur->is_node());
#endif
      // Find all the local sub-entities to be send to other rank
      std::vector<hcell_t *> cells;
      traversal(
        cur,
        [&](hcell_t * cell, std::vector<hcell_t *> & c) {
          if(cell->key() == cur->key())
            return true;
          if(cell->is_node()) {
            c.push_back(cell);
            return true;
          }
          c.push_back(cell);
          return false;
        } // lambda
        ,
        cells);

      for(int j = 0; j < cells.size(); ++j) {
        if(cells[j]->is_node()) {
          cells[j]->set_nchildren_to_receive(cells[j]->nchildren());
          tmp_nodes_replies.emplace_back(
            cells[j]->owner(), cells[j]->key(), *get_node(cells[j]),
            cells[j]->nchildren());
        }
        else if(cells[j]->is_entity()) {
          tmp_entities_replies.emplace_back(
            cells[j]->owner(), cells[j]->key(), *get_entity(cells[j]));
        }
#ifdef _DEBUG_TREE_
        else {
          assert(false);
        } // if
#endif
      } // for
    } // for
    if(tmp_nodes_replies.size() != 0) {
      mpi_replies_[current_replies_].push_back(MPI_Request{});
      nodes_replies_.push_back(tmp_nodes_replies);
      MPI_Issend(&nodes_replies_[nodes_replies_.size() - 1][0],
        sizeof(share_node_t) * tmp_nodes_replies.size(), MPI_BYTE, partner,
        REPLY_NODE, MPI_COMM_WORLD, &mpi_replies_[current_replies_].back());
      found = true;
      if(mpi_replies_[current_replies_].size() >= requests_keys_max_ - 1) {
        current_replies_++;
        mpi_replies_.resize(current_replies_ + 1);
        mpi_replies_[current_replies_].reserve(requests_keys_max_);
      } // if
    } // if
    if(tmp_entities_replies.size() != 0) {
      mpi_replies_[current_replies_].push_back(MPI_Request{});
      entities_replies_.push_back(tmp_entities_replies);
      MPI_Issend(&entities_replies_[entities_replies_.size() - 1][0],
        sizeof(share_entity_t) * tmp_entities_replies.size(), MPI_BYTE, partner,
        REPLY_ENTITY, MPI_COMM_WORLD, &mpi_replies_[current_replies_].back());
      found = true;
      if(mpi_replies_[current_replies_].size() >= requests_keys_max_ - 1) {
        current_replies_++;
        mpi_replies_.resize(current_replies_ + 1);
        mpi_replies_[current_replies_].reserve(requests_keys_max_);
      } // if
    } // if
#ifdef _DEBUG_TREE_
    assert(found);
#endif
  }

  /**
   * @brief Receive a group of entities from a distant node and
   * add them in the local tree, create the appropriate parents
   * of the entity inserted to be able to reach it
   */
  void recv_entity_replies_(const int & partner, const int & nrecv) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nentities = nrecv / sizeof(share_entity_t);
    std::vector<share_entity_t> recv_entities(nentities);
    MPI_Recv(&recv_entities[0], nrecv, MPI_BYTE, partner, REPLY_ENTITY,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(int i = 0; i < recv_entities.size(); ++i) {
      key_t pkey = recv_entities[i].key;
      pkey.pop();
      auto parent = htable_.find(pkey);
      shared_entities_.push_back(recv_entities[i].entity);
#ifdef _DEBUG_TREE_
      assert(htable_.find(recv_entities[i].key) == htable_.end());
#endif
      htable_.emplace(recv_entities[i].key,
        hcell_t(recv_entities[i].key, shared_entities_.size() - 1));
      auto it = htable_.find(recv_entities[i].key);
      it->second.set_shared();
      it->second.set_owner(recv_entities[i].owner);
      // Change parent
      int child = recv_entities[i].key.last_value();
      parent->second.add_child(child);
    } // for
  }

  /**
   * @brief Add a group of nodes received from another rank in the
   * local tree. This may request the creation of intermediate nodes
   * in the tree.
   */
  void recv_node_replies_(const int & partner, const int & nrecv) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nnodes = nrecv / sizeof(share_node_t);
    std::vector<share_node_t> recv_nodes(nnodes);
    MPI_Recv(&recv_nodes[0], nrecv, MPI_BYTE, partner, REPLY_NODE,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for(int i = 0; i < nnodes; ++i) {
      key_t pkey = recv_nodes[i].key;
      pkey.pop();
      auto parent = htable_.find(pkey);
#ifdef _DEBUG_TREE_
      assert(parent != htable_.end());
#endif
      shared_nodes_.push_back(recv_nodes[i].node);
#ifdef _DEBUG_TREE_
      assert(htable_.find(recv_nodes[i].key) == htable_.end());
#endif
      htable_.emplace(recv_nodes[i].key, recv_nodes[i].key);
      auto it = htable_.find(recv_nodes[i].key);
      it->second.set_shared();
      it->second.set_node_idx(shared_nodes_.size() - 1);
      it->second.set_owner(recv_nodes[i].owner);
      it->second.set_nchildren_to_receive(recv_nodes[i].nchildren);
      // Change parent
      int child = recv_nodes[i].key.last_value();
#ifdef _DEBUG_TREE_
      key_t ckey = recv_nodes[i].key;
      ckey.pop();
      assert(parent->first == ckey);
#endif
      parent->second.add_child(child);
    } // for
    // Do we need to clean a node after being requested
    // since it should never be requested again.
    // Otherwise we need to store:
    // current children + total children (bitmap)
    // if(parent->second.nchildren() ==
    // get_node(&parent->second)->sub_entities())
    //  parent->second.unset_requested();
  }

  /**
   * @brief Share the entities/nodes with neighbors
   * Find the branches that are not allocated yet, hey are on the limit
   * of the domain.
   */
  template<typename CCOFM>
  void share_nodes_(CCOFM && f_cc) {
    double start = omp_get_wtime();
    log_one(trace) << "Sharing nodes/entities " << std::endl;

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;
    // Do the hypercube communciation to share the branches
    // Add them in the tree in the same time
    int dim = log2(size);
    bool non_power_2 = false;
    // In case of non power two, consider dim + 1
    // In this case all rank also take size + 1 rank
    int ghosts_rank = -1;
    if(1 << dim < size) {
      non_power_2 = true;
      dim++;
      ghosts_rank = rank + (1 << (dim - 1));
      if(ghosts_rank > (1 << dim) - 1)
        ghosts_rank = rank;
    }
    std::vector<share_entity_t> ghosts_entities, r_ghosts_entities;
    std::vector<share_node_t> ghosts_nodes, r_ghosts_nodes;
    const int sz_entities = sizeof(share_entity_t);
    const int sz_nodes = sizeof(share_node_t);
    int s_ge_size, s_gn_size;

    // For non power of 2 cases, store the array of entities
    // of the ghost rank represented.
    std::vector<share_entity_t> r_ghosts_entities_n2;
    std::vector<share_node_t> r_ghosts_nodes_n2;

    // Communication for all channels
    for(int i = 0; i < dim; ++i) {
      int partner = rank ^ (1 << i);
      if(partner >= 0 && partner < size) {
        ghosts_entities.clear();
        ghosts_nodes.clear();
        find_nodes_(ghosts_nodes, ghosts_entities);
#ifdef _DEBUG_TREE_
        assert(partner >= 0 && partner != rank && partner < size);
#endif
        // Send lobound and hibound and bytes for nodes/entities
        s_ge_size = ghosts_entities.size() * sz_entities;
        s_gn_size = ghosts_nodes.size() * sz_nodes;
        std::pair<int[2], key_t[2]> s_keys;
        s_keys.first[0] = s_ge_size;
        s_keys.first[1] = s_gn_size;
        s_keys.second[0] = lobound_;
        s_keys.second[1] = hibound_;
        std::pair<int[2], key_t[2]> s_rkeys;
        MPI_Sendrecv(&s_keys, sizeof(std::pair<int, key_t[2]>), MPI_BYTE,
          partner, 0, &s_rkeys, sizeof(std::pair<int, key_t[2]>), MPI_BYTE,
          partner, 0, MPI_COMM_WORLD, &status);
        lobound_ = std::min(s_rkeys.second[0], lobound_);
        hibound_ = std::max(s_rkeys.second[1], hibound_);
        // Send entities
        r_ghosts_entities.resize(s_rkeys.first[0] / sz_entities);
        MPI_Sendrecv(&ghosts_entities[0], s_ge_size, MPI_BYTE, partner, 0,
          &r_ghosts_entities[0], s_rkeys.first[0], MPI_BYTE, partner, 0,
          MPI_COMM_WORLD, &status);
        // Send nodes
        r_ghosts_nodes.resize(s_rkeys.first[1] / sz_nodes);
        MPI_Sendrecv(&ghosts_nodes[0], s_gn_size, MPI_BYTE, partner, 0,
          &r_ghosts_nodes[0], s_rkeys.first[1], MPI_BYTE, partner, 0,
          MPI_COMM_WORLD, &status);
        // Insert the nodes/entities in the tree
        for(size_t j = 0; j < r_ghosts_entities.size(); ++j) {
          if(r_ghosts_entities[j].owner != rank) {
            shared_entities_.push_back(r_ghosts_entities[j].entity);
            load_shared_entity_(shared_entities_.size() - 1,
              r_ghosts_entities[j].key, r_ghosts_entities[j].owner);
          }
        }
        for(size_t j = 0; j < r_ghosts_nodes.size(); ++j) {
          if(r_ghosts_nodes[j].owner != rank) {
            shared_nodes_.push_back(r_ghosts_nodes[j].node);
            load_shared_node_(shared_nodes_.size() - 1, r_ghosts_nodes[j].key,
              r_ghosts_nodes[j].owner);
          }
        }
        cofm_update_(root(), f_cc);
        // Handle the non power two cases
      } // if
      if(non_power_2) {
        // If this ghosts_rank exists for a real rank don't use it
        if(rank != ghosts_rank && ghosts_rank <= size - 1)
          continue;
        int ghosts_partner = ghosts_rank ^ (1 << i);
        partner = ghosts_partner;
        // Case already handled before
        if(ghosts_rank == rank && partner < size)
          continue;
        if(ghosts_rank > size && partner > size)
          continue;
        if(partner >= size) {
          partner -= (1 << (dim - 1));
        }
        if(partner == rank) {
          // Add into each buffer, no communication needed
          ghosts_entities.clear();
          ghosts_nodes.clear();
          find_nodes_(ghosts_nodes, ghosts_entities);
          r_ghosts_entities_n2.insert(r_ghosts_entities_n2.end(),
            ghosts_entities.begin(), ghosts_entities.end());
          r_ghosts_nodes_n2.insert(
            r_ghosts_nodes_n2.end(), ghosts_nodes.begin(), ghosts_nodes.end());
          for(size_t j = 0; j < r_ghosts_entities_n2.size(); ++j) {
            if(r_ghosts_entities_n2[j].owner != rank) {
              shared_entities_.push_back(r_ghosts_entities_n2[j].entity);
              load_shared_entity_(shared_entities_.size() - 1,
                r_ghosts_entities_n2[j].key, r_ghosts_entities_n2[j].owner);
            }
          }
          for(size_t j = 0; j < r_ghosts_nodes_n2.size(); ++j) {
            if(r_ghosts_nodes_n2[j].owner != rank) {
              shared_nodes_.push_back(r_ghosts_nodes_n2[j].node);
              load_shared_node_(shared_nodes_.size() - 1,
                r_ghosts_nodes_n2[j].key, r_ghosts_nodes_n2[j].owner);
            }
          }
          cofm_update_(root(), f_cc);
        }
        else {
          if(ghosts_rank == rank) {
            ghosts_entities.clear();
            ghosts_nodes.clear();
            find_nodes_(ghosts_nodes, ghosts_entities);
          }
          else {
            ghosts_entities = r_ghosts_entities_n2;
            ghosts_nodes = r_ghosts_nodes_n2;
          }
#ifdef _DEBUG_TREE_
          assert(partner >= 0 && partner != rank && partner < size);
#endif
          // Send lobound and hibound and bytes for nodes/entities
          s_ge_size = ghosts_entities.size() * sz_entities;
          s_gn_size = ghosts_nodes.size() * sz_nodes;
          std::pair<int[2], key_t[2]> s_keys;
          s_keys.first[0] = s_ge_size;
          s_keys.first[1] = s_gn_size;
          s_keys.second[0] = lobound_;
          s_keys.second[1] = hibound_;
          std::pair<int[2], key_t[2]> s_rkeys;
          MPI_Sendrecv(&s_keys, sizeof(std::pair<int, key_t[2]>), MPI_BYTE,
            partner, 0, &s_rkeys, sizeof(std::pair<int, key_t[2]>), MPI_BYTE,
            partner, 0, MPI_COMM_WORLD, &status);
          lobound_ = std::min(s_rkeys.second[0], lobound_);
          hibound_ = std::max(s_rkeys.second[1], hibound_);
          // Send entities
          r_ghosts_entities.resize(s_rkeys.first[0] / sz_entities);
          MPI_Sendrecv(&ghosts_entities[0], s_ge_size, MPI_BYTE, partner, 0,
            &r_ghosts_entities[0], s_rkeys.first[0], MPI_BYTE, partner, 0,
            MPI_COMM_WORLD, &status);
          // Send nodes
          r_ghosts_nodes.resize(s_rkeys.first[1] / sz_nodes);
          MPI_Sendrecv(&ghosts_nodes[0], s_gn_size, MPI_BYTE, partner, 0,
            &r_ghosts_nodes[0], s_rkeys.first[1], MPI_BYTE, partner, 0,
            MPI_COMM_WORLD, &status);
          if(rank == ghosts_rank) {
            // Insert the nodes/entities in the tree
            for(size_t j = 0; j < r_ghosts_entities.size(); ++j) {
              if(r_ghosts_entities[j].owner != rank) {
                shared_entities_.push_back(r_ghosts_entities[j].entity);
                load_shared_entity_(shared_entities_.size() - 1,
                  r_ghosts_entities[j].key, r_ghosts_entities[j].owner);
              }
            }
            for(size_t j = 0; j < r_ghosts_nodes.size(); ++j) {
              if(r_ghosts_nodes[j].owner != rank) {
                shared_nodes_.push_back(r_ghosts_nodes[j].node);
                load_shared_node_(shared_nodes_.size() - 1,
                  r_ghosts_nodes[j].key, r_ghosts_nodes[j].owner);
              }
            }
            cofm_update_(root(), f_cc);
          }
          else {
            r_ghosts_entities_n2.insert(r_ghosts_entities_n2.end(),
              ghosts_entities.begin(), ghosts_entities.end());
            r_ghosts_nodes_n2.insert(r_ghosts_nodes_n2.end(),
              ghosts_nodes.begin(), ghosts_nodes.end());
          }
        } // if
      } // if
    } // for
#ifdef _DEBUG_TREE_
    assert(root()->is_node());
#endif
    log_one(trace) << "Sharing nodes/entities.done: " << omp_get_wtime() - start
                   << "s" << std::endl;
  }

  /**
   * @brief Complete the CoFM in the tree with new entities and branches
   * This function is called during the sharing of entities/nodes.
   * The inserted parents in the tree needs to be associated with a
   * cofm and it needs to be computed.
   */
  template<typename CCOFM>
  void cofm_update_(hcell_t * current, CCOFM && f_cc) {
    key_t nkey = current->key();
    key_t min_key, max_key;
    key_boundary_(nkey, min_key, max_key);
    if(min_key >= lobound_ && max_key <= hibound_) {
      if(current->is_unset()) {
        std::vector<hcell_t *> daughters;
        daughters.reserve(nchildren_);
        for(int i = 0; i < nchildren_; ++i) {
          if(current->get_child(i)) {
            key_t ckey = nkey;
            ckey.push(i);
            auto it = htable_.end(); 
            it = htable_.find(ckey);
#ifdef _DEBUG_TREE_
            assert(it != htable_.end());
#endif
            daughters.push_back(&(htable_.find(ckey)->second));
          } // if
        } // for
        for(size_t i = 0; i < daughters.size(); ++i) {
          cofm_update_(daughters[i], f_cc);
        } // for
        current->set_shared();
        current->set_node_idx(shared_nodes_.size());
        shared_nodes_.push_back(nkey);
        cofm_children_(&shared_nodes_[current->node_idx()], daughters, f_cc);
      } // if
    }
  }

  /**
   * @brief Find the nodes or entities to be shared with other ranks.
   * This searches for the first nodes/entities that have an index.
   */
  void find_nodes_(std::vector<share_node_t> & nodes,
    std::vector<share_entity_t> & entities) {
    nodes.clear();
    entities.clear();
    std::vector<hcell_t *> queue;
    std::vector<hcell_t *> nqueue;
    queue.push_back(root());
    while(!queue.empty()) {
      for(hcell_t * cur : queue) {
        key_t nkey = cur->key();
        if(cur->is_unset()) {
#ifdef _DEBUG_TREE_
          assert(cur->type() != 0);
#endif
          for(int j = 0; j < nchildren_; ++j) {
            if(cur->get_child(j)) {
              key_t ckey = nkey;
              ckey.push(j);
              auto it = htable_.find(ckey);
              nqueue.push_back(&(it->second));
            } // if
          } // for
        }
        else {
          if(cur->is_node()) {
            cofm_t * cofm = get_node(cur);
            // TODO: check if initializing nchildren with 0 is OK here
            nodes.emplace_back(cur->owner(), cur->key(), *cofm, 0); 
          }
          else {
            entity_t * ent = get_entity(cur);
            entities.emplace_back(cur->owner(), cur->key(), *ent);
          } // if
        } // else
      } // for
      queue.clear();
      queue = nqueue;
      nqueue.clear();
    } // while
  }

  /**
   * @brief Load an entity in the tree from a distant process
   * Call the add_parent_ function to link this entity to
   * the tree.
   **/
  void
  load_shared_entity_(const int & entity_idx, key_t key, const int & owner) {
#ifdef _DEBUG_TREE_
    assert(htable_.find(key) == htable_.end());
#endif
    htable_.emplace(key, hcell_t(key, entity_idx));
    hcell_t * cur = &(htable_.find(key)->second);
    cur->set_shared();
    cur->set_owner(owner);
    int lastbit = key.pop_value();
    add_parent_(key, lastbit, owner);
  }

  /**
   * @brief Load a node in the tree from a distant process
   * Call the add_parent_ function to link this entity to
   * the tree
   **/
  void load_shared_node_(const int & node_idx, key_t key, const int & owner) {
#ifdef _DEBUG_TREE_
    assert(htable_.find(key) == htable_.end());
#endif
    htable_.emplace(key, key);
    hcell_t * cur = &(htable_.find(key)->second);
    cur->set_shared();
    cur->set_node_idx(node_idx);
    cur->set_owner(owner);
    int lastbit = key.pop_value();
    add_parent_(key, lastbit, owner);
  }

  /**
   * @brief Add missing parent in the tree from a distant
   * node or entity insertion.
   * The new parents are empty and we add the children in
   * the types for the tree traversal
   */
  void add_parent_(key_t key, int child, const int & owner) {
    auto parent = htable_.end();
    while((parent = htable_.find(key)) == htable_.end()) {
      htable_.emplace(key, key);
      parent = htable_.find(key);
      parent->second.set_shared();
      parent->second.add_child(child);
      parent->second.set_owner(owner);
      child = key.pop_value();
    } // while
#ifdef _DEBUG_TREE_
    assert(parent->second.node_idx() == -1);
#endif
    parent->second.add_child(child);
  }

  /**
   * @brief Finish a branch during the creation of the tree
   * This branch is done and this rank is the only one that
   * have information for it (middle of its local tree).
   * Add node cofm + compute cofm data with cofm_children_
   */
  template<typename CCOFM>
  void finish_(const key_t & key, CCOFM && f_c) {
    hcell_t * n = &(htable_.find(key)->second);
    if(n->entity_idx() != -1)
      return;
#ifdef _DEBUG_TREE_
    assert(n->node_idx() == -1 && n->entity_idx() == -1);
#endif
    n->set_node_idx(cofm_.size());
    cofm_.emplace_back(key);
    std::vector<hcell_t *> daughters;
    daughters.reserve(nchildren_);
    for(int j = 0; j < nchildren_; ++j) {
      if(n->get_child(j)) {
        key_t ckey = key;
        ckey.push(j);
        auto it = htable_.end(); 
        it = htable_.find(ckey);
#ifdef _DEBUG_TREE_
        assert(it != htable_.end());
#endif
        daughters.push_back(&(htable_.find(ckey)->second));
      } // if
    } // for
    cofm_children_(&cofm_[n->node_idx()], daughters, f_c);
  }

  /**
   * @brief Return a pointer to the hcell daughters of a node.
   * Using the key of the current hcell and pushing the child number in
   * the parent key.
   */
  void daughters_(hcell_t * cell, hcell_t ** daughters, int & children) {
#ifdef _DEBUG_TREE_
    assert(cell != nullptr);
    assert(daughters != nullptr);
#endif
    key_t nkey = cell->key();
    children = 0;
    for(int i = 0; i < nchildren_; ++i) {
      if(cell->get_child(i)) {
        key_t ckey = nkey;
        ckey.push(i);
        auto it = htable_.find(ckey);
#ifdef _DEBUG_TREE_
        assert(it != htable_.end());
#endif
        daughters[children++] = &(it->second);
      } // if
    } // for
  }

  /**
   * @brief Compute the CofM data based on the daughters of the node.
   */
  template<typename CCOFM>
  void cofm_children_(cofm_t * cofm,
    const std::vector<hcell_t *> & daughters,
    CCOFM && f_ce) {
    std::vector<entity_t *> v_entities;
    std::vector<cofm_t *> v_nodes;
    for(size_t i = 0; i < daughters.size(); ++i) {
      if(daughters[i]->is_entity()) {
        v_entities.push_back(get_entity(daughters[i]));
      }
      else if(daughters[i]->is_node()) {
        v_nodes.push_back(get_node(daughters[i]));
      }
#ifdef _DEBUG_TREE_
      else {
        assert(false);
      }
#endif
    } // for
    // Compute center of mass values
    f_ce(cofm, v_entities, v_nodes);
  }

  /**
   * @brief Exchange the keys boundries with direct neighbors
   * This will be used during the sharing of nodes and entities.
   * The first and last process have: 100.. and 177.. for low or high.
   **/
  void exchange_boundaries_(const key_t & hikey,
    const key_t & lokey,
    key_t & hibound,
    key_t & lobound) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    int up = rank == size - 1 ? -1 : rank + 1;
    int down = rank == 0 ? -1 : rank - 1;

    if(rank & 1) {
      if(up >= 0) {
        MPI_Sendrecv(&hikey, sizeof(key_t), MPI_BYTE, up, 0, &hibound,
          sizeof(key_t), MPI_BYTE, up, 0, MPI_COMM_WORLD, &status);
      }
      else {
        hibound = key_t::max();
      }
      if(down >= 0) {
        MPI_Sendrecv(&lokey, sizeof(key_t), MPI_BYTE, down, 0, &lobound,
          sizeof(key_t), MPI_BYTE, down, 0, MPI_COMM_WORLD, &status);
      }
      else {
        lobound = key_t::min();
      }
    }
    else {
      if(up >= 0) {
        MPI_Sendrecv(&hikey, sizeof(key_t), MPI_BYTE, up, 0, &hibound,
          sizeof(key_t), MPI_BYTE, up, 0, MPI_COMM_WORLD, &status);
      }
      else {
        hibound = key_t::max();
      }
      if(down >= 0) {
        MPI_Sendrecv(&lokey, sizeof(key_t), MPI_BYTE, down, 0, &lobound,
          sizeof(key_t), MPI_BYTE, down, 0, MPI_COMM_WORLD, &status);
      }
      else {
        lobound = key_t::min();
      }
    } // if
  } // exchange_boundaries_

  /**
   * @brief Find the min (1XX00..) and max (1XX77..) keys around a key
   * in the tree. This key might not be at the max depth.
   */
  void key_boundary_(const key_t & key, key_t & min_key, key_t & max_key) {
    key_t stop;
    stop = key_t::min();
    max_key = key;
    min_key = key;
    while(min_key < stop) {
      max_key.push((1 << dimension) - 1);
      min_key.push(0);
    } // while
  } // key_boundary

  /**
   * @brief Initialization of the communication array
   * based on the number of ranks involved in the comms.
   */
  void init_comms_(const int & size) {
    std::fill(comms_done_.begin(), comms_done_.end(), false);
    mpi_requests_.resize(1);
    mpi_requests_[0].reserve(requests_keys_max_);
    mpi_replies_.resize(1);
    mpi_replies_[0].reserve(requests_keys_max_);
    current_requests_ = 0;
    current_replies_ = 0;
    comms_timer_ = 0;
    lost_timer_ = 0;
  }

  /**
   * @brief Clean the communications arrays.
   * Check the requests termination
   */
  void clean_comms_() {
    // Check that all communications have beens completed
    MPI_Status status;
    int flag;
    for(int i = 0; i < mpi_requests_.size(); ++i) {
      for(int j = 0; j < mpi_requests_[i].size(); ++j) {
        MPI_Test(&mpi_requests_[i][j], &flag, &status);
#ifdef _DEBUG_TREE_
        assert(flag);
#endif
      }
    }
    mpi_requests_.clear();
    for(int i = 0; i < mpi_replies_.size(); ++i) {
      for(int j = 0; j < mpi_replies_[i].size(); ++j) {
        MPI_Test(&mpi_replies_[i][j], &flag, &status);
#ifdef _DEBUG_TREE_
        assert(flag);
#endif
      }
    }
    mpi_replies_.clear();
    requests_keys_.clear();
    nodes_replies_.clear();
    entities_replies_.clear();
  }

  // KEEP this hashing function to be able to
  // switch from hashtable to unordered_map from
  // std.
  template<class key_t>
  struct branch_id_hasher__ {
    size_t operator()(const key_t & k) const noexcept {
      return static_cast<size_t>(k.value() & ((1 << 22) - 1));
    }
  };

  // Tree topology
  size_t max_depth_;
  // KEEP this to switch with hashtable
  // to see the best implementation
  using umap_t = std::unordered_map<key_t, hcell_t, branch_id_hasher__<key_t>>;
  // using umap_t = hashtable<key_t, hcell_t>;
  typename umap_t::iterator root_;
  umap_t htable_;
  range_t range_;
  std::vector<cofm_t> cofm_;
  std::vector<entity_t> entities_;
  std::vector<entity_t> shared_entities_;
  std::vector<cofm_t> shared_nodes_;
  static constexpr int nchildren_ = (1 << dimension);
  key_t hibound_, lobound_;
  // Communication
  std::vector<std::vector<key_t>> requests_keys_;
  int current_requests_, current_replies_;
  std::vector<std::vector<MPI_Request>> mpi_requests_;
  std::vector<std::vector<MPI_Request>> mpi_replies_;
  std::vector<std::vector<share_node_t>> nodes_replies_;
  std::vector<std::vector<share_entity_t>> entities_replies_;
  std::vector<bool> comms_done_;
  bool comms_all_done_;
  const int requests_keys_max_ = 100;
  double comms_timer_, lost_timer_;
  // Traversal
  const int sub_entities_ = 128;
  const int fmm_sub_entities_ = 0;
};

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
