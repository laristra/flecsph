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

#include "tree.h"
#include "utils.h"

#include "params.h" // For the variable smoothing length 

using namespace mpi_utils;

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

  // To share the ghosts data within the radius
  mpi_ghosts_t ghosts_data;

  const int criterion_branches = 64; // Number of branches to share in 
                                       // exchange_branches
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
    // double start = omp_get_wtime(); // unused

    // Sort the keys 
    std::sort(rbodies.begin(),rbodies.end(),
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
    
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<std::pair<entity_key_t,int64_t>> splitters;
    
    generate_splitters_samples(splitters,rbodies,totalnbodies);

    // The rbodies array is already sorted. We just need to determine the 
    // limits for each process
    // Reset the class buffer 
    reset_buffers();
    
    int cur_proc = 0;

    for(auto bi: rbodies){
      if(cur_proc == size-1){
        // Last process reached, just copy 
        scount[cur_proc]++; 
        continue; 
      } // if
      if(bi.first < splitters[cur_proc].first){
        scount[cur_proc]++;
      }else{
        scount[++cur_proc]++;
      } // if
    } // for

    std::vector<std::pair<entity_key_t,body>> recvbuffer; 
    mpi_alltoallv(scount,rbodies,recvbuffer); 
    
    rbodies.clear();
    rbodies = recvbuffer; 
    
    // Sort the incoming buffer 
    sort(rbodies.begin(),rbodies.end(),
      [](auto& left, auto &right){
        if(left.first<right.first){
          return true; 
        }
        if(left.first == right.first){
          return left.second.getId()<right.second.getId(); 
        }
        return false; 
      }); // sort
   

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
      std::vector<body_holder_mpi_t> tmpsendbuffer; 
      for(int j = cur; j < cur+count[i]; ++j ){
        // Then for each branches
        tree_topology_t::subentity_space_t ents; 
        if(param::sph_variable_h){
          ents = tree.find_in_box(recv_branches[j][0],recv_branches[j][1],
            tree_geometry_t::intersects_sphere_box);
        }else{
          ents = tree.find_in_box(recv_branches[j][0],recv_branches[j][1],
            tree_geometry_t::within_box);
        } 
        for(auto ent: ents){
          assert(ent != nullptr);
          tmpsendbuffer.push_back(body_holder_mpi_t{
            ent->coordinates(),rank,ent->mass(),ent->getBody()->getId(),
            ent->getBody()->getSmoothinglength()});
        }
      }

      std::sort(tmpsendbuffer.begin(),tmpsendbuffer.end(),
        [](const auto& left, const auto& right){
          return left.id < right.id;
      });  

      tmpsendbuffer.erase(
          std::unique(tmpsendbuffer.begin(),tmpsendbuffer.end(),
          [](const auto& left, const auto& right){
            return (left.id == right.id);
          }),tmpsendbuffer.end());
      scount[i] = tmpsendbuffer.size(); 

      sendbuffer.insert(sendbuffer.end(),tmpsendbuffer.begin(),
        tmpsendbuffer.end());

      cur += count[i];
    }

    std::vector<body_holder_mpi_t> recvbuffer;
    mpi_alltoallv(scount,sendbuffer,recvbuffer); 

    // Add them in the tree 
    for(auto bi: recvbuffer)
    {
      assert(bi.owner!=rank);
      assert(bi.mass!=0.);
      auto nbi = tree.make_entity(bi.position,nullptr,bi.owner,bi.mass,bi.id,
          bi.h);
      tree.insert(nbi);
      assert(nbi->global_id() == bi.id);
    }

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
    auto itsb = ghosts_data.sbodies.begin();
    for(auto bi: ghosts_data.sholders)
    {
      assert(bi->getBody()!=nullptr);
      *itsb = *(bi->getBody());
      itsb++;
    }  

    MPI_Alltoallv(&ghosts_data.sbodies[0],&ghosts_data.nsbodies[0],
      &ghosts_data.soffsets[0],MPI_BYTE,
      &ghosts_data.rbodies[0],&ghosts_data.nrbodies[0],
      &ghosts_data.roffsets[0],MPI_BYTE,MPI_COMM_WORLD);

    // Then link the holders with these bodies
    auto ents = tree.entities().to_vec(); 
    int64_t nelem = ents.size(); 
    
    int64_t totalfound = 0;
#pragma omp parallel for reduction(+:totalfound)
    for(int64_t i=0; i<nelem; ++i)
    {
      body_holder* bi = ents[i];
      if(!bi->is_local())
      {
        for(auto& nl: ghosts_data.rbodies)
        {
          if(bi->global_id() == nl.id())
          {
            totalfound++;
            bi->setBody(&(nl));
            break;
          }
        }
      }
    }

    if(totalfound != (int) ghosts_data.rbodies.size())
      std::cout<<rank<<": t="<<totalfound<<" g="<<ghosts_data.rbodies.size()<<
       std::endl; 
    
    assert(totalfound == (int)ghosts_data.rbodies.size()); 
  
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
    tree_topology_t& tree,
    std::vector<body_holder*>& lbodies
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

    int64_t nelem = lbodies.size();

    // Count send
#pragma omp parallel for 
    for(int64_t i=0; i<nelem;++i)
    {
      // array of bools to check unique send 
      std::vector<bool> proc(size,false);
      proc[rank] = true; 
      body_holder * bi = lbodies[i];
      assert(bi->is_local());
      tree_topology_t::subentity_space_t nbs; 
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
      std::vector<bool> proc(size,false);
      proc[rank] = true;
      body_holder* bi = lbodies[i];
      assert(bi->is_local());
      tree_topology_t::subentity_space_t nbs; 
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

    int64_t totalsendbodies = 0L; 
    int64_t totalrecvbodies = 0L; 
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
    for(int i=0;i<size;++i){
      ghosts_data.nsbodies[i]*=sizeof(body);
      ghosts_data.nrbodies[i]*=sizeof(body);
      ghosts_data.soffsets[i]*=sizeof(body);
      ghosts_data.roffsets[i]*=sizeof(body);
    }
  
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
          
    for(auto bi: bodies){
      for(size_t i=0;i<gdimension;++i){
        if(bi.second.coordinates()[i]+bi.second.getSmoothinglength()>
            lrange[1][i])
          lrange[1][i] = bi.second.coordinates()[i]+
                        bi.second.getSmoothinglength();
        if(bi.second.coordinates()[i]-bi.second.getSmoothinglength()<
            lrange[0][i])
          lrange[0][i] = bi.second.coordinates()[i]-
                        bi.second.getSmoothinglength();
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
    size_t maxnsamples = 1024;
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

    // Generate the splitters
    splitters.resize(size-1);
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

      //std::cout<<entity_key_t::first_key()<<std::endl;
      for(int i=0;i<size-1;++i){
        int64_t position = (master_nkeys/size)*(i+1);
        splitters[i] = master_keys[position];
        //std::cout<<splitters[i].first<<std::endl;
      } // for
      // Last key 
      //std::cout<<entity_key_t::last_key()<<std::endl;
    } // if

    // Bradcast the splitters 
    MPI_Bcast(&splitters[0],(size-1)*sizeof(std::pair<entity_key_t,int64_t>)
    ,MPI_BYTE,0,MPI_COMM_WORLD);
  }

}; // class tree_colorer 

#endif // _mpisph_tree_colorer_h_
