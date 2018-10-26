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
 * @file physics.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic physics implementation
 */

#ifndef _boundary_h_
#define _boundary_h_

#include <vector>

#include "params.h"

namespace boundary{
  using namespace param;
 
  point_t min_boundary; 
  point_t max_boundary;
  double damp; 


  void pboundary_clean(
      std::vector<std::pair<entity_key_t,body>>& lbodies)
  {
    // Delete all local WALL particles
    for(auto it = lbodies.begin(); it != lbodies.end(); )
    {
      if((*it).second.getType() == particle_type_t::WALL)
        it = lbodies.erase(it);
      else
        ++it;
    }
  }

  void pboundary_generate(
      std::vector<std::pair<entity_key_t,body>>& lbodies,
      double halo_size)
  {
    int64_t original_nparticles = lbodies.size();
    MPI_Allreduce(MPI_IN_PLACE,&original_nparticles,1,MPI_INT64_T,
        MPI_SUM,MPI_COMM_WORLD); 

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
        point_t coord = lbodies[i].second.coordinates();
        if(lbodies[i].second.coordinates()[d] > box[1][d])
        {
          coord[d] = box[0][d] + (coord[d]-box[1][d]); 
        }else if(lbodies[i].second.coordinates()[d] < box[0][d])
        {
          coord[d] = box[1][d] - (box[0][d]-coord[d]);
        }
        lbodies[i].second.setPosition(coord);
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
        if(lbodies[i].second.coordinates()[d]+halo_size> box[1][d] ||
          lbodies[i].second.coordinates()[d]-halo_size< box[0][d])
        {
          on_edge[d] = true; 

          // New on this axis 
          body nu = lbodies[i].second; 
          nu.setType(particle_type_t::WALL); 
          point_t coord = nu.coordinates(); 
          if(nu.coordinates()[d]+halo_size>box[1][d]){
            coord[d] = box[0][d] - fabs((nu.coordinates()[d]+box[0][d]));
          }else{
            coord[d] = box[1][d] + fabs((nu.coordinates()[d]+box[1][d]));
          }
          nu.setPosition(coord); 
          #pragma omp critical 
            edge.push_back(nu);
        }
      }
      // If several dimensions interfere, add in the corners 
      if(gdimension == 1)
        continue; 
      if(gdimension > 1){
        if(on_edge[0] && on_edge[1]){
          body nu = lbodies[i].second; 
          nu.setType(particle_type_t::WALL);
          point_t coord = nu.coordinates(); 
          if(nu.coordinates()[0]+halo_size>box[1][0]){
            coord[0] = box[0][0] - fabs(nu.coordinates()[0]+box[0][0]);
          }else{
            coord[0] = box[1][0] + fabs(nu.coordinates()[0]+box[1][0]);
          }
           if(nu.coordinates()[1]+halo_size>box[1][1]){
            coord[1] = box[0][1] - fabs(nu.coordinates()[1]+box[0][1]);
          }else{
            coord[1] = box[1][1] + fabs(nu.coordinates()[1]+box[1][1]);
          }
          nu.setPosition(coord);
          #pragma omp critical 
            edge.push_back(nu);
        }
      }
    }
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    //std::cout<<rank<<" "<<edge.size()<<" particles transfered"<<
    //std::endl<<std::flush; 
    
    // Add them in the body array
    for(auto b: edge)
      lbodies.push_back(std::make_pair(entity_key_t{},b));
  
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
      if(b.second.getType() == WALL)
        b.second.setId((i++));
      total_check += b.second.id(); 
    }
    MPI_Allreduce(MPI_IN_PLACE,&total_check,1,MPI_INT64_T,MPI_SUM,
        MPI_COMM_WORLD);
    //std::cout<<total_check<<" == "<<(total_part-1)*(total_part)/2<<std::endl;
    assert(total_check == (total_part-1)*(total_part)/2);
  }
    
  /**
   * @brief      Apply boundaries if they are set
   *
   * @param      srch  The source's body holder
   *
   * @return     True if the particle have been considered outside the 
   * boundaries
   */
  bool
  compute_boundaries(
      body_holder* srch)
  {
    body* source = srch->getBody();
    point_t velocity = source->getVelocity();
    point_t position = source->getPosition();
    point_t velocityHalf = source->getVelocityhalf();

    bool considered = false;

    if(stop_boundaries){
      bool stop = false; 
      for(size_t i = 0; i < gdimension; ++i){
        if(position[i] < min_boundary[i] ||
          position[i] > max_boundary[i]){
          stop = true; 
        }
      }
      if(stop){
        velocity = point_t{};
        velocityHalf = point_t{};
        considered = true;
      
      }
    }else if(reflect_boundaries){
      for(size_t dim=0;dim < gdimension ; ++dim){
        if(position[dim] < min_boundary[dim] || 
            position[dim] > max_boundary[dim]){
          double barrier = max_boundary[dim];
          if(position[dim] < min_boundary[dim]){
            barrier = min_boundary[dim];
          }

          // Here just invert the velocity vector and velocityHalf 
          double tbounce = (position[dim]-barrier)/velocity[dim];
          position -= velocity*(1-damp)*tbounce;

          position[dim] = 2*barrier-position[dim];
          velocity[dim] = -velocity[dim];
          velocityHalf[dim] = -velocityHalf[dim];

          velocity *= damp;
          velocityHalf *= damp;
          considered = true;
        }
      }
    }
    source->setPosition(position);
    source->setVelocity(velocity);
    source->setVelocityhalf(velocityHalf);
    return considered;
  }

}; // physics

#endif // _default_physics_h_
