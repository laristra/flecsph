/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file physics.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic physics implementation
 */

#ifndef _boundary_h_
#define _boundary_h_

#include <vector>

#ifdef BOOST
#include <boost/sort/sort.hpp>
//#include <boost/compute/algorithm/stable_partition.hpp>
#endif 

#include "params.h"

namespace boundary{
  using namespace param;

  point_t min_boundary;
  point_t max_boundary;
  double damp;


  void pboundary_clean(
      std::vector<body>& lbodies)
  {
    // Sort items based on type, all the walls at the end
    auto start = std::stable_partition(
        lbodies.begin(),lbodies.end(),
        [](const auto& left)
        {
          return left.getType() != particle_type_t::WALL;
        });

    assert(start != lbodies.begin());

#ifdef DEBUG
    for(auto it = start; it != lbodies.end(); ++it)
      assert(it->getType() == particle_type_t::WALL);
#endif

    lbodies.erase(start,lbodies.end());

    // Delete all the last ones

    // Delete all local WALL particles
    //for(auto it = lbodies.begin(); it != lbodies.end(); )
    //{
    //  if((*it).second.getType() == particle_type_t::WALL)
    //    it = lbodies.erase(it);
    //  else
    //    ++it;
    //}
  }

  void
  add_corner_XX(
    const std::array<bool,gdimension>& on_edge,
    const std::array<size_t,2>& d,
    const std::array<point_t,2>& box,
    const double& halo_size,
    const body& lbody,
    std::vector<body>& edge
  )
  {
    if(on_edge[d[0]] && on_edge[d[1]]){
      body nu = lbody;
      nu.setType(particle_type_t::WALL);
      point_t coord = nu.coordinates();
      for(size_t i = 0 ; i < 2 ; ++i){
        if(nu.coordinates()[d[i]]+halo_size>box[1][d[i]]){
          coord[d[i]] = box[0][d[i]] - abs(nu.coordinates()[d[i]]+box[0][d[i]]);
        }else{
          coord[d[i]] = box[1][d[i]] + abs(nu.coordinates()[d[i]]+box[1][d[i]]);
        }
      }
      nu.set_coordinates(coord);
      #pragma omp critical
        edge.push_back(nu);
    }
  }

  void
  add_corner_XXX(
    const std::array<bool,gdimension>& on_edge,
    const std::array<size_t,3>& d,
    const std::array<point_t,2>& box,
    const double& halo_size,
    const body& lbody,
    std::vector<body>& edge
  )
  {
    if(on_edge[d[0]] && on_edge[d[1]] && on_edge[d[2]]){
      body nu = lbody;
      nu.setType(particle_type_t::WALL);
      point_t coord = nu.coordinates();
      for(size_t i = 0 ; i < 3 ; ++i){
        if(nu.coordinates()[d[i]]+halo_size>box[1][d[i]]){
          coord[d[i]] = box[0][d[i]] - abs(nu.coordinates()[d[i]]+box[0][d[i]]);
        }else{
          coord[d[i]] = box[1][d[i]] + abs(nu.coordinates()[d[i]]+box[1][d[i]]);
        }
      }
      nu.set_coordinates(coord);
      #pragma omp critical
        edge.push_back(nu);
    }
  }

  void pboundary_generate(
      std::vector<body>& lbodies,
      double halo_size)
  {
    int64_t original_nparticles = lbodies.size();
    MPI_Allreduce(MPI_IN_PLACE,&original_nparticles,1,MPI_INT64_T,
        MPI_SUM,MPI_COMM_WORLD);

    bool periodic[3] = {param::periodic_boundary_x,
      param::periodic_boundary_y,
      param::periodic_boundary_z};

    // Generate box from dimension
    range_t box;
    point_t min, max;
    // At least 1D
    min[0] = -box_length/2.; max[0] = box_length/2.;
    if(gdimension >= 1){
      min[1] = -box_width/2.; max[1] = box_width/2.;
    }
    if(gdimension >= 2){
      min[2] = -box_height/2.; max[2] = box_height/2.;
    }
    box = {min,max};


    // Step 1, teleport the particles to the other side of the domain
    // Keep the same id
    #pragma omp parallel for
    for(int i = 0 ; i < lbodies.size(); ++i)
    {
      for(int d = 0 ; d < gdimension ; ++d)
      {
        // If not in this axis, do not move the particles
        if(!periodic[d]) continue;

        point_t coord = lbodies[i].coordinates();
        if(lbodies[i].coordinates()[d] > box[1][d])
        {
          coord[d] = box[0][d] + (coord[d]-box[1][d]);
        }else if(lbodies[i].coordinates()[d] < box[0][d])
        {
          coord[d] = box[1][d] - (box[0][d]-coord[d]);
        }
        lbodies[i].set_coordinates(coord);
      }
    }

    // Search for particles near the edge
    std::vector<body> edge;
    #pragma omp parallel for
    for(int i = 0; i < lbodies.size(); ++i)
    {
      std::array<bool,gdimension> on_edge;
      on_edge.fill(false);

      for(int d = 0; d < gdimension ; ++d)
      {
        if(!periodic[d]) continue;

        if(lbodies[i].coordinates()[d]+halo_size> box[1][d] ||
          lbodies[i].coordinates()[d]-halo_size< box[0][d])
        {
          on_edge[d] = true;

          // New on this axis
          body nu = lbodies[i];
          nu.setType(particle_type_t::WALL);
          point_t coord = nu.coordinates();
          if(nu.coordinates()[d]+halo_size>box[1][d]){
            coord[d] = box[0][d] - fabs((nu.coordinates()[d]+box[0][d]));
          }else{
            coord[d] = box[1][d] + fabs((nu.coordinates()[d]+box[1][d]));
          }
          nu.set_coordinates(coord);
          #pragma omp critical
            edge.push_back(nu);
        }
      }
      // If several dimensions interfere, add in the corners
      if(gdimension == 1)
        continue;
      if(gdimension == 2){
        // x and y
        add_corner_XX(on_edge,{0,1},box,halo_size,lbodies[i],edge);
      }
      if(gdimension == 3){
        // Multiple combinations
        // xy, xz, yz, xyz
        add_corner_XX(on_edge,{0,1},box,halo_size,lbodies[i],edge);
        add_corner_XX(on_edge,{0,2},box,halo_size,lbodies[i],edge);
        add_corner_XX(on_edge,{2,1},box,halo_size,lbodies[i],edge);

        add_corner_XXX(on_edge,{0,1,2},box,halo_size,lbodies[i],edge);
      }
    }
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    //std::cout<<rank<<" "<<edge.size()<<" particles transfered"<<
    //std::endl<<std::flush;

    // Add them in the body array
    for(auto b: edge)
      lbodies.push_back(b);

    std::vector<int64_t> nparticles_wall(size);
    int64_t local_wall = edge.size();
    MPI_Allgather(&local_wall,1,MPI_INT64_T,
        &nparticles_wall[0],1,MPI_INT64_T,MPI_COMM_WORLD);

    // Prefix scan
    std::partial_sum(nparticles_wall.begin(),
        nparticles_wall.end(),nparticles_wall.begin());
    nparticles_wall.insert(nparticles_wall.begin(),0);

    // Re Index the particles in this case to avoid conflicts
    // Keep the same index for domain particles
    // Create special index for wall particles
    int64_t i = original_nparticles+nparticles_wall[rank];
    int64_t total_check = 0L;
    int64_t total_part = original_nparticles+nparticles_wall[size];
    for(auto& b: lbodies)
    {
      if(b.getType() == WALL)
        b.set_id((i++));
      total_check += b.id();
    }
    MPI_Allreduce(MPI_IN_PLACE,&total_check,1,MPI_INT64_T,MPI_SUM,
        MPI_COMM_WORLD);
    //std::cout<<total_check<<" == "<<(total_part-1)*(total_part)/2<<std::endl;
    assert(total_check == (total_part-1)*(total_part)/2);
  }

}; // physics

#endif // _default_physics_h_
