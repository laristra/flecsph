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
 * @file tree_colorer.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Function needed for MPI distribution of the bodies
 */

#ifndef _mpisph_tree_colorer_h_
#define _mpisph_tree_colorer_h_

#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <omp.h>

#include <boost/sort/sort.hpp>

#include "tree.h"
#include "utils.h"
#include "default_physics.h"

#include "params.h" // For the variable smoothing length

using namespace mpi_utils;

//#define BOOST_PARALLEL 1
// Output the data regarding the distribution for debug
#define OUTPUT_TREE_INFO 1


/**
 * @brief       Structure to keep the data during the ghosts sharing.
 * Fill the structure during compute_ghosts and then exchange during
 * refresh_ghosts ''
 */
struct mpi_ghosts_t{
  std::vector<body> sbodies;
  std::vector<body> rbodies;
  std::vector<int>  nsbodies;
  std::vector<int>  soffsets;
  std::vector<int>  nrbodies;
  std::vector<int>  roffsets;
  std::vector<body_holder*> sholders;
  std::vector<body_holder*> rholders;
};

struct body_holder_mpi_t{
  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::point__<element_t, dimension>;

  point_t position;
  int owner;
  double mass;
  flecsi::topology::entity_id_t id;
  double h;
  entity_key_t key;
};

/**
 * @brief      All the function and buffers for the tree_colorer.
 *
 * @tparam     T     Type of the class
 * @tparam     D     Dimension for the problem
 * @todo fix it for type
 */
template<
  typename T,
  size_t D
>
class tree_colorer
{
private:
  // For all the MPI_Alltoallv communications
  std::vector<int> rcount;
  std::vector<int> scount;
  std::vector<int> roffset;
  std::vector<int> soffset;

  // Communication of body_holders
  MPI_Datatype MPI_BH_T;

  // To share the ghosts data within the radius
  mpi_ghosts_t ghosts_data;

  const int criterion_branches = 1; // Number of sub-entities in the branches
  const size_t noct = 256*1024;        // Number of octets used for quicksort

  void reset_buffers()
  {
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    rcount.resize(size);
    scount.resize(size);
    roffset.resize(size);
    soffset.resize(size);

    ghosts_data.sbodies.resize(size);
    ghosts_data.rbodies.resize(size);
    ghosts_data.nsbodies.resize(size);
    ghosts_data.nrbodies.resize(size);
    ghosts_data.roffsets.resize(size);
    ghosts_data.soffsets.resize(size);

    std::fill(rcount.begin(),rcount.end(),0);
    std::fill(scount.begin(),scount.end(),0);
    std::fill(roffset.begin(),roffset.end(),0);
    std::fill(soffset.begin(),soffset.end(),0);
  }

public:
  static const size_t dimension = D;
  using point_t = flecsi::point__<T,dimension>;


  tree_colorer(){
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    rcount.resize(size);
    scount.resize(size);
    roffset.resize(size);
    soffset.resize(size);

    ghosts_data.sbodies.resize(size);
    ghosts_data.rbodies.resize(size);
    ghosts_data.nsbodies.resize(size);
    ghosts_data.nrbodies.resize(size);
    ghosts_data.roffsets.resize(size);
    ghosts_data.soffsets.resize(size);

    // Create the types for communications
    // Create datatype
    body_holder_mpi_t tmp;
    MPI_Datatype type[5] =
      {MPI_DOUBLE,MPI_INT,MPI_DOUBLE,MPI_INT64_T,MPI_DOUBLE};
    int blocklen[5] = { gdimension, 2, 1, 1, 2};
    MPI_Aint disp[5];
    disp[0] = 0;
    disp[1] = sizeof(double)*gdimension;
    disp[2] = disp[1] + sizeof(int)*2; // Times 2 due to alignement
    disp[3] = disp[2] + sizeof(double);
    disp[4] = disp[3] + sizeof(int64_t);
    MPI_Type_create_struct(5, blocklen, disp, type, &MPI_BH_T);
    MPI_Type_commit(&MPI_BH_T);

    int size_BH;
    MPI_Type_size(MPI_BH_T,&size_BH);
    rank|| clog(trace) << "Size of MPI_BH_T: "<<size_BH<< " BH: "<<sizeof(body_holder_mpi_t)<<std::endl;
  }

  ~tree_colorer(){}

/*~---------------------------------------------------------------------------*
 * Function for sorting and distribution
 *~---------------------------------------------------------------------------*/

  /**
  * @brief      Sorting of the input particles or current particles using MPI.
  * This method is composed of several steps to implement the quick sort:
  * - Each process sorts its local particles
  * - Each process generate a subset of particles to fit the byte size limit
  * to send
  * - Each subset if send to the master (Here 0) who generates the pivot for
  * quick sort
  * - Pivot are send to each processes and they create buckets based on the
  * pivot
  * - Each bucket is send to the owner
  * - Each process sorts its local particles again
  * This is the first implementation, it can be long for the first sorting but
  * but then as the particles does not move very fast, the particles on the
  * edge are the only ones shared
  *
  * @param      rbodies       The rbodies, local bodies of this process.
  * @param[in]  totalnbodies  The totalnbodies on the overall simulation.
  */
  void mpi_qsort(
    std::vector<std::pair<entity_key_t,body>>& rbodies,
    int totalnbodies)
  {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Sort the keys
    // Use boost parallel sort
#ifdef BOOST_PARALLEL
    boost::sort::block_indirect_sort(
#else
    std::sort(
#endif
        rbodies.begin(),rbodies.end(),
        [](auto& left, auto& right){
        if(left.first < right.first){
          return true;
        }
        if(left.first == right.first){
          return left.second.getId() < right.second.getId();
        }
        return false;
      }); // sort

    // If one process, done
    if(size==1){
      clog(trace)<<"Local particles: "<<totalnbodies<<std::endl;
      return;
    } // if

    std::vector<std::pair<entity_key_t,int64_t>> splitters;
    generate_splitters_samples(splitters,rbodies,totalnbodies);

    // The rbodies array is already sorted. We just need to determine the
    // limits for each process
    // Reset the class buffer
    reset_buffers();

    int cur_proc = 0;

    assert(splitters.size() == size-1+2);

    int64_t nbodies = rbodies.size();
    for(size_t i = 0L ; i < nbodies; ++i){
      if(rbodies[i].first >= splitters[cur_proc].first &&
        rbodies[i].first < splitters[cur_proc+1].first){
          scount[cur_proc]++;
        }else{
          i--;
          cur_proc++;
        }
    }

    // Check that we considered all the bodies
    assert( std::accumulate(scount.begin(), scount.end(), 0) == rbodies.size());


    //for(auto bi: rbodies){
    //  if(cur_proc >= size-1){
        // Last process reached, just copy
    //    scount[cur_proc]++;
    //    continue;
    //  } // if
    //  if(bi.first < splitters[cur_proc].first){
    //    scount[cur_proc]++;
    //  }else{
    //    while(!(bi.first < splitters[cur_proc].first)){
    //      cur_proc++;
    //    }
    //    scount[cur_proc]++;
    //  } // if
    //} // for

    std::vector<std::pair<entity_key_t,body>> recvbuffer;

    // Direct exchange using point to point
    mpi_alltoallv_p2p(scount,rbodies,recvbuffer);

    rbodies.clear();
    rbodies = recvbuffer;

#ifdef OUTPUT
    std::vector<int> totalprocbodies;
    totalprocbodies.resize(size);
    int mybodies = rbodies.size();
    // Share the final array size of everybody
    MPI_Allgather(&mybodies,1,MPI_INT,&totalprocbodies[0],1,MPI_INT,
      MPI_COMM_WORLD);
    #ifdef OUTPUT_TREE_INFO
    std::ostringstream oss;
    oss<<"Repartition: ";
    for(auto num: totalprocbodies)
      oss<<num<<";";
    double end = omp_get_wtime();
    rank|| clog(trace)<<oss.str()<<std::endl;
    #endif
#endif // OUTPUT
  } // mpi_qsort

  /**
   * @brief      Exchange the useful branches of the current tree of the procs.
   * There is several ways to share. Here we look for the particles in the
   * global bounding box and share them. The branches are constructed directly
   * by the FleCSI tree structure. A better way should be to share branches
   * instead. But we need to change the tree structure in FleCSI for that.
   *
   * This function provide the non local particles that are use to find the
   * ghosts.
   *
   * @param      tree             The tree
   * @param      rbodies          The rbodies, local bodies of this process
   * @param      ranges           The ranges, the key range of all the processes
   * @param      range_total      The range total of all the particles
   */
  void
  mpi_branches_exchange(
    tree_topology_t& tree,
    std::vector<std::pair<entity_key_t,body>>& rbodies,
    std::vector<std::array<point_t,2>>& ranges,
    std::array<point_t,2>& range_total
  )
  {
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef OUTPUT
    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef OUTPUT_TREE_INFO
    rank|| clog(trace)<<"Branches repartition" << std::endl << std::flush;
    #endif
#endif

    // Do a tree search up to a branch
    // Keep those branches in a list
    std::vector<branch_t*> search_branches;
    tree.find_sub_cells(
      tree.root(),
      criterion_branches,
      search_branches);

    rank|| clog(trace) << "1. branches: "<< search_branches.size() << std::endl << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);

    // Make a list of boundaries
    std::vector<range_t> send_branches(search_branches.size());

    #pragma omp parallel for
    for(long unsigned int i=0;i<search_branches.size();++i){
      send_branches[i][0] = search_branches[i]->bmin();
      send_branches[i][1] = search_branches[i]->bmax();
    }

    std::vector<range_t> recv_branches;
    std::vector<int> count;
    mpi_allgatherv(send_branches,recv_branches,count);

    // Output branches
    //if(rank == 0)
    //  output_branches_VTK(recv_branches,count,physics::iteration);

    rank|| clog(trace) << "2. Received:"<<recv_branches.size() << std::endl << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);


    reset_buffers();
    std::vector<body_holder_mpi_t> sendbuffer;
    int cur = 0;

    // Compute the requested nodes
    // Search in the tree for each processes
    for(int i=0;i<size;++i)
    {
      if(i==rank){
        cur += count[i];
        continue;
      }

      // Predicate for the sets
      auto set_predicate = [](const auto& left, const auto& right)->bool
      {
        return left.id < right.id;
      };

      std::set<body_holder_mpi_t,decltype(set_predicate)> tmpsendbuffer(set_predicate);

      #pragma omp parallel shared(tmpsendbuffer)
      {
        // Local set for the threads
        std::set<body_holder_mpi_t,decltype(set_predicate)> tmpsendset(set_predicate);

        #pragma omp for
        for(int j = cur; j < cur+count[i]; ++j ){
          // Then for each branches
          tree_topology_t::entity_space_ptr_t ents;
          if(param::sph_variable_h){
            ents = tree.find_in_box(recv_branches[j][0],recv_branches[j][1],
              tree_geometry_t::intersects_sphere_box);
          }else{
            ents = tree.find_in_box(recv_branches[j][0],recv_branches[j][1],
              tree_geometry_t::within_box);
          }

          for(auto ent: ents){
            // Mark these bodies as shared for the future
            ent->set_shared();
            assert(ent != nullptr);
            tmpsendset.insert(body_holder_mpi_t{
              ent->coordinates(),rank,ent->mass(),ent->getBody()->getId(),
              ent->getBody()->getSmoothinglength(),ent->get_entity_key()});
          }
        } // for
        // Merge the sets
        #pragma omp critical
          tmpsendbuffer.merge(tmpsendset);
      } // pragma omp parallel

      scount[i] = tmpsendbuffer.size();

      sendbuffer.insert(sendbuffer.end(),tmpsendbuffer.begin(),
        tmpsendbuffer.end());

      cur += count[i];
    }

    rank|| clog(trace) << "3. Entities:"<< sendbuffer.size() << std::endl<<std::flush;
    rank|| clog(trace) << "BH:"<<sizeof(body_holder_mpi_t)<<" tot:"<<sendbuffer.size()<<std::endl<<std::flush;
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<body_holder_mpi_t> recvbuffer;

    std::vector<int> recvcount(size), recvoffsets(size), sendoffsets(size);
    // Exchange the send count
    MPI_Alltoall(&scount[0],1,MPI_INT,&recvcount[0],1,MPI_INT,
        MPI_COMM_WORLD);
    std::partial_sum(recvcount.begin(),recvcount.end(),&recvoffsets[0]);
    recvoffsets.insert(recvoffsets.begin(),0);
    std::partial_sum(scount.begin(),scount.end(),&sendoffsets[0]);
    sendoffsets.insert(sendoffsets.begin(),0);
    // Set the recvbuffer to the right size
    recvbuffer.resize(recvoffsets.back());
    // Transform the offsets for bytes
    #pragma omp parallel for
    for(int i=0;i<size;++i){
      assert(scount[i]>=0);
      assert(recvcount[i]>=0);
      assert(sendoffsets[i]>=0);
      assert(recvoffsets[i]>=0);
    } // for
    std::vector<MPI_Status> status(size);
    std::vector<MPI_Request> request(size);

#pragma omp parallel
{
    #pragma omp for nowait
    for(int i = 0 ; i < size; ++i){
      if(scount[i] != 0){
        MPI_Isend(&(sendbuffer[sendoffsets[i]]),scount[i],MPI_BH_T,
          i,0,MPI_COMM_WORLD,&request[i]);
      }
    }
    #pragma omp for nowait
    for(int i = 0 ; i < size; ++i){
      if(recvcount[i] != 0){
        MPI_Recv(&(recvbuffer[recvoffsets[i]]),recvcount[i],MPI_BH_T,
          i,MPI_ANY_TAG,MPI_COMM_WORLD,&status[i]);
        // Add in the tree
        #pragma omp critical
        {
          for(size_t j = recvoffsets[i]; j < recvoffsets[i]+recvcount[i]; ++j )
          {
            auto* bi = &(recvbuffer[j]);
            assert(bi->owner!=rank);
            assert(bi->mass!=0.);
            auto id = tree.make_entity(bi->key,bi->position,nullptr,bi->owner,
              bi->mass,bi->id,bi->h);
            tree.insert(id);
            auto nbi = tree.get(id);
            assert(!nbi->is_local());
            assert(nbi->global_id() == bi->id);
          } // for
        } // omp critical
      }
      if(scount[i] != 0){
        MPI_Wait(&request[i],&status[i]);
      }
    }
} // omp parallel

    rank|| clog(trace) << "4. Ent received:"<<recvbuffer.size()<<std::endl<<std::flush;
    //MPI_Barrier(MPI_COMM_WORLD);

    // Add them in the tree
    // Not doable in parallel due to the tree utilization
    //for(auto bi: recvbuffer)
    //{
    //  assert(bi.owner!=rank);
    //  assert(bi.mass!=0.);
    //  auto id = tree.make_entity(bi.position,nullptr,bi.owner,bi.mass,bi.id,
    //      bi.h);
    //  tree.insert(id);
    //  auto nbi = tree.get(id);
    //  assert(!nbi->is_local());
    //  assert(nbi->global_id() == bi.id);
    //}



#ifdef OUTPUT_TREE_INFO
    MPI_Barrier(MPI_COMM_WORLD);
    rank || clog(trace)<<".done "<<std::endl;
#endif

  }


void mpi_refresh_ghosts(
    tree_topology_t& tree
    )
  {
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if(size == 1){
      return;
    }

#ifdef OUTPUT_TREE_INFO
    MPI_Barrier(MPI_COMM_WORLD);
    rank|| clog(trace)<<"Refresh Ghosts" << std::flush;
    double start = omp_get_wtime();
#endif
   // Refresh the sendbodies with new data
   // auto itsb = ghosts_data.sbodies.begin();

    #pragma omp parallel for
    for(size_t i = 0; i < ghosts_data.sbodies.size() ; ++i){
      assert(ghosts_data.sholders[i]->getBody() != nullptr);
      ghosts_data.sbodies[i] = *(ghosts_data.sholders[i]->getBody());
    }

    // Point to point communication
    std::vector<MPI_Status> status(size);
    std::vector<MPI_Request> request(size);
#pragma omp parallel
{
    #pragma omp for nowait
    for(int i = 0 ; i < size; ++i){
      if(ghosts_data.nsbodies[i] != 0){
        char * start = (char*)&(ghosts_data.sbodies[0]);
        MPI_Isend(start+ghosts_data.soffsets[i],ghosts_data.nsbodies[i],
            MPI_BYTE,i,0,MPI_COMM_WORLD,&request[i]);
      }
    } // for
    #pragma omp for nowait
    for(int i = 0 ; i < size; ++i){
      if(ghosts_data.nrbodies[i] != 0){
        char * start = (char*)&(ghosts_data.rbodies[0]);
        MPI_Recv(start+ghosts_data.roffsets[i],ghosts_data.nrbodies[i],
            MPI_BYTE,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status[i]);
      }
      // Link the received ghosts
      #pragma omp parallel for
      for(size_t j = ghosts_data.roffsets[i]/sizeof(body); j <
         ghosts_data.roffsets[i]/sizeof(body)+
         ghosts_data.nrbodies[i]/sizeof(body) ; ++j){
        auto * bh = tree.get_ghost(ghosts_data.rbodies[j].id());
        assert(!bh->is_local());
        bh->setBody(&(ghosts_data.rbodies[j]));
      } // for
      if(ghosts_data.nsbodies[i] != 0){
        MPI_Wait(&request[i],&status[i]);
      }
    } // for
} // omp parallel

    //MPI_Alltoallv(&ghosts_data.sbodies[0],&ghosts_data.nsbodies[0],
    //  &ghosts_data.soffsets[0],MPI_BYTE,
    //  &ghosts_data.rbodies[0],&ghosts_data.nrbodies[0],
    //  &ghosts_data.roffsets[0],MPI_BYTE,MPI_COMM_WORLD);

//#pragma omp parallel for
  //For all the received neighbors, need to find a local ghosts
  //for(size_t i = 0 ; i < ghosts_data.rbodies.size(); ++i){
  //  auto * bh = tree.get_ghost(ghosts_data.rbodies[i].id());
  //  assert(!bh->is_local());
  //  bh->setBody(&(ghosts_data.rbodies[i]));
  //}

#ifdef OUTPUT_TREE_INFO
    MPI_Barrier(MPI_COMM_WORLD);
    rank|| clog(trace) <<".done "<< std::endl << std::flush;

#endif
  }


/**
   * @brief      Prepare the buffer for the ghost transfer function.
   * Based on the non local particles shared in the mpi_branches_exchange,
   * this function extract the really needed particles and find the ghosts.
   * Then, as those ghosts can be requested several times in an iteration,
   * the buffer are set and can bne use in mpi_refresh_ghosts.
   *
   * @param      tree             The tree
   * @param      range            The range
   */
  void
  mpi_compute_ghosts(
    tree_topology_t& tree
  )
  {
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // No need to compute ghosts for one process
    if(size == 1){
      return;
    }

#ifdef OUTPUT_TREE_INFO
    MPI_Barrier(MPI_COMM_WORLD);
    rank|| clog(trace)<<"Compute Ghosts" << std::flush;
    double start = omp_get_wtime();
#endif

    // Clean the structure
    ghosts_data.sbodies.clear();
    ghosts_data.rbodies.clear();
    ghosts_data.nsbodies.clear();
    ghosts_data.nrbodies.clear();

    ghosts_data.nsbodies.resize(size);
    ghosts_data.nrbodies.resize(size);
    std::fill(ghosts_data.nsbodies.begin(),ghosts_data.nsbodies.end(),0);
    std::fill(ghosts_data.nrbodies.begin(),ghosts_data.nrbodies.end(),0);

    int64_t nelem = tree.entities().size();

    // Count send
#pragma omp parallel for
    for(int64_t i=0; i<nelem;++i)
    {

      body_holder * bi = tree.get(i);
      if(!bi->is_shared()) continue;

      // array of bools to check unique send
      std::vector<bool> proc(size,false);
      proc[rank] = true;

      assert(bi->is_local());
      tree_topology_t::entity_space_ptr_t nbs;
      if(param::sph_variable_h){
        nbs = tree.find_in_radius(
            bi->coordinates(),
            bi->getBody()->getSmoothinglength(),
            tree_geometry_t::within_square
        );
      }else{
        nbs = tree.find_in_radius(
            bi->coordinates(),
            bi->getBody()->getSmoothinglength(),
            tree_geometry_t::within
        );
      }
      for(auto nb: nbs)
      {
        if(!nb->is_local() && !proc[nb->owner()])
        {
          // Mark this particle as sent for this process
          proc[nb->owner()] = true;
#pragma omp atomic update
          ghosts_data.nsbodies[nb->owner()]++;
        } // if
      } // for
    } // for

    int64_t totalsbodies=0;
    // Total
    for(int i=0;i<size;++i)
    {
      totalsbodies += ghosts_data.nsbodies[i];
    }

    std::vector<int> offset(size,0);
    for(int i=1; i<size; ++i)
    {
      offset[i] += offset[i-1]+ghosts_data.nsbodies[i-1];
    }

    assert(totalsbodies>=0);
    // Allocate the send array
    ghosts_data.sbodies.resize(totalsbodies);
    ghosts_data.sholders.resize(totalsbodies);
    // Temp variable to offset in the sbodies array
    std::vector<int> spbodies(size,0);
    // Fill the vector
#pragma omp parallel for
    for(int64_t i=0; i<nelem; ++i)
    {

      body_holder* bi = tree.get(i);
      if(!bi->is_local()) continue;

      std::vector<bool> proc(size,false);
      proc[rank] = true;

      assert(bi->is_local());
      tree_topology_t::entity_space_ptr_t nbs;
      if(param::sph_variable_h){
        nbs = tree.find_in_radius(
            bi->coordinates(),
            bi->getBody()->getSmoothinglength(),
            tree_geometry_t::within_square
        );
      }else{
        nbs = tree.find_in_radius(
            bi->coordinates(),
            bi->getBody()->getSmoothinglength(),
            tree_geometry_t::within
        );
      }
      for(auto nb: nbs)
      {
        if(!nb->is_local() && !proc[nb->owner()])
        {
          proc[nb->owner()] = true;
          int pos = 0;

#pragma omp atomic capture
          pos = spbodies[nb->owner()]++;

          // Write
          pos += offset[nb->owner()];
          assert(pos<totalsbodies);
          ghosts_data.sholders[pos] = bi;
          ghosts_data.sbodies[pos] = *(bi->getBody());
        } // if
      } // for
    } // for

    MPI_Alltoall(&ghosts_data.nsbodies[0],1,MPI_INT,
        &ghosts_data.nrbodies[0],1,MPI_INT,MPI_COMM_WORLD);

#ifdef DEBUG
    // total receive
    int64_t totalnrecv = 0;
    for(size_t i = 0 ; i < ghosts_data.nrbodies.size(); ++i)
      totalnrecv += ghosts_data.nrbodies[i];
    std::vector<int64_t> tabnrecv(size);
    MPI_Gather(&totalnrecv,1,MPI_INT64_T,&(tabnrecv[0]),1,MPI_INT64_T,0,
        MPI_COMM_WORLD);
    if(rank == 0){
      std::ostringstream oss;
      for(size_t i = 0 ; i < size; ++i){
        oss << tabnrecv[i] << ";";
      }
      clog(trace) << oss.str() << std::endl;
    }
#endif

    int64_t totalsendbodies = 0L;
    int64_t totalrecvbodies = 0L;

#pragma omp parallel for reduction(+:totalsendbodies) \
    reduction(+:totalrecvbodies)
    for(int i=0;i<size;++i){
      assert(ghosts_data.nsbodies[i]>=0);
      assert(ghosts_data.nrbodies[i]>=0);
      totalsendbodies += ghosts_data.nsbodies[i];
      totalrecvbodies += ghosts_data.nrbodies[i];
    } // for

    // Prepare offsets for alltoallv
    ghosts_data.roffsets[0]=0;
    ghosts_data.soffsets[0]=0;

    for(int i=1;i<size;++i){
      ghosts_data.roffsets[i] = ghosts_data.nrbodies[i-1]+
        ghosts_data.roffsets[i-1];
      ghosts_data.soffsets[i] = ghosts_data.nsbodies[i-1]+
        ghosts_data.soffsets[i-1];
    }

    ghosts_data.rbodies.resize(totalrecvbodies);

    // Convert the offsets to byte
#pragma omp parallel for
    for(int i=0;i<size;++i){
      ghosts_data.nsbodies[i]*=sizeof(body);
      ghosts_data.nrbodies[i]*=sizeof(body);
      ghosts_data.soffsets[i]*=sizeof(body);
      ghosts_data.roffsets[i]*=sizeof(body);
    }

    assert(totalsendbodies == ghosts_data.sbodies.size());

#ifdef OUTPUT_TREE_INFO
    MPI_Barrier(MPI_COMM_WORLD);
    double end = omp_get_wtime();
    rank|| clog(trace)<<".done "<< end-start << "s"<<std::endl;
#endif
  }

/*~---------------------------------------------------------------------------*
 * Utils functions
 *~---------------------------------------------------------------------------*/

  /**
   * @brief      Compute the global range of all the particle system
   *
   * @param      bodies           The bodies, local of this process
   * @param      range            The range computed in that function
   */
  void
  mpi_compute_range(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    std::array<point_t,2>& range)
  {
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Compute the local range
    range_t lrange;

    lrange[1] = bodies.back().second.coordinates();
    lrange[0] = bodies.back().second.coordinates();

#pragma omp parallel
{
  range_t trange;
  trange[1] = bodies.back().second.coordinates();
  trange[0] = bodies.back().second.coordinates();

#pragma omp parallel for
    for(size_t i = 0 ; i < bodies.size(); ++i){
      for(size_t d = 0; d < gdimension; ++d){

        if(bodies[i].second.coordinates()[d]+
            bodies[i].second.getSmoothinglength()>trange[1][d])
          trange[1][d] = bodies[i].second.coordinates()[d]+
                        bodies[i].second.getSmoothinglength();

        if(bodies[i].second.coordinates()[d]-
            bodies[i].second.getSmoothinglength()<trange[0][d])
          trange[0][d] = bodies[i].second.coordinates()[d]-
                        bodies[i].second.getSmoothinglength();
      }
    }
#pragma omp critical
    for(size_t d = 0 ; d < gdimension; ++d)
    {
      lrange[1][d] = std::max(lrange[1][d],trange[1][d]);
      lrange[0][d] = std::min(lrange[0][d],trange[0][d]);
    }
}

    double max[gdimension];
    double min[gdimension];
    for(size_t i=0; i<gdimension; ++i){
      max[i] = lrange[1][i];
      min[i] = lrange[0][i];
    }

    // Do the MPI Reduction
    MPI_Allreduce(MPI_IN_PLACE,max,dimension,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,min,dimension,MPI_DOUBLE,MPI_MIN,
        MPI_COMM_WORLD);

    for(size_t d = 0; d < gdimension ; ++d){
      range[0][d] = min[d];
      range[1][d] = max[d];
    }
  }


/**
 * @brief      Use in mpi_qsort to generate the splitters to sort the particles
 * in the quick sort algorithm
 * In this function we take some samplers of the total particles and the root
 * determines the splitters
 * This version is based on the sample splitter algorithm but we generate more
 * samples on each process
 *
 * @param      splitters  The splitters used in the qsort in mpi_qsort
 * @param[in]  rbodies  The local bodies of the process
 */
  void
  generate_splitters_samples(
    std::vector<std::pair<entity_key_t,int64_t>>& splitters,
    std::vector<std::pair<entity_key_t,body>>& rbodies,
    const int64_t totalnbodies
  ){
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Create a vector for the samplers
    std::vector<std::pair<entity_key_t,int64_t>> keys_sample;
    // Number of elements for sampling
    // In this implementation we share up to 256KB to
    // the master.
    size_t maxnsamples = noct/sizeof(std::pair<entity_key_t,int64_t>);
    int64_t nvalues = rbodies.size();
    size_t nsample = maxnsamples*((double)nvalues/(double)totalnbodies);

    if(nvalues<(int64_t)nsample){nsample = nvalues;}

    for(size_t i=0;i<nsample;++i){
      int64_t position = (nvalues/(nsample+1.))*(i+1.);
      keys_sample.push_back(std::make_pair(rbodies[position].first,
      rbodies[position].second.getId()));
    } // for
    assert(keys_sample.size()==(size_t)nsample);

    std::vector<std::pair<entity_key_t,int64_t>> master_keys;
    std::vector<int> master_recvcounts;
    std::vector<int> master_offsets;
    int master_nkeys = 0;

    if(rank==0){
      master_recvcounts.resize(size);
    } // if

    // Echange the number of samples
    MPI_Gather(&nsample,1,MPI_INT,
      &master_recvcounts[0],1,MPI_INT,0,MPI_COMM_WORLD);

    // Master
    // Sort the received keys and create the pivots
    if(rank == 0){
      master_offsets.resize(size);
      master_nkeys = std::accumulate(
        master_recvcounts.begin(),master_recvcounts.end(),0);
      if(totalnbodies<master_nkeys){master_nkeys=totalnbodies;}
      // Number to receiv from each process
      for(int i=0;i<size;++i){
        master_recvcounts[i]*=sizeof(std::pair<entity_key_t,int64_t>);
      } // for
      std::partial_sum(master_recvcounts.begin(),master_recvcounts.end(),
        &master_offsets[0]);
      master_offsets.insert(master_offsets.begin(),0);
      master_keys.resize(master_nkeys);
    } // if

    MPI_Gatherv(&keys_sample[0],nsample*sizeof(std::pair<entity_key_t,int64_t>)
      ,MPI_BYTE,&master_keys[0],&master_recvcounts[0],&master_offsets[0]
      ,MPI_BYTE,0,MPI_COMM_WORLD);

    // Generate the splitters, add zero and max keys
    splitters.resize(size-1+2);
    if(rank==0){
      std::sort(master_keys.begin(),master_keys.end(),
        [](auto& left, auto& right){
          if(left.first < right.first){
            return true;
          }
          if(left.first == right.first){
            return left.second < right.second;
          }
          return false;
        });

      splitters[0].first = entity_key_t::min();
      splitters[0].second = 0L;
      splitters[size].first = entity_key_t::max();
      splitters[size].second = LONG_MAX;

      for(int i=0;i<size-1;++i){
        int64_t position = (master_nkeys/size)*(i+1);
        splitters[i+1] = master_keys[position];
        assert(splitters[i+1].first > splitters[0].first &&
          splitters[i+1].first < splitters[size].first);
      } // for

      // Print the keys
      //for(auto k: splitters){
      //  std::cout<<k.first<<std::endl;
      //}
    } // if

    // Bradcast the splitters
    MPI_Bcast(&splitters[0],(size-1+2)*sizeof(std::pair<entity_key_t,int64_t>)
    ,MPI_BYTE,0,MPI_COMM_WORLD);
  }

}; // class tree_colorer

#endif // _mpisph_tree_colorer_h_
