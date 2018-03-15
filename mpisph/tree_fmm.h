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
 * @file tree_fmm.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Function needed for FMM computation
 */

#ifndef _mpisph_tree_fmm_h_
#define _mpisph_tree_fmm_h_

#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <omp.h>

#include "tree.h"

template<
  typename T,
  size_t D
>
class tree_fmm
{

  //using point_t = flecsi::point<T,D>;
  using geometry_t = flecsi::topology::tree_geometry<T,D>;
  static const size_t dimension = D;

  /**
   * @brief      Structure for COM communication during the gravitation
   *             computation process
  */
  struct mpi_cell_t{
    point_t position;
    point_t fc = {};
    double dfcdr[dimension*dimension] = {};
    double dfcdrdr[dimension*dimension*dimension] = {};
    point_t bmax;
    point_t bmin;
    branch_id_t id;
    int ninterations = 0;
  
    mpi_cell_t(
      point_t position_,
      point_t bmin_,
      point_t bmax_,
      branch_id_t id_)
    {
      position = position_;
      id = id_;
      bmax = bmax_;
      bmin = bmin_;
    };

    mpi_cell_t(){};
  };

private:

  // Used in the COM computation and exchange
  std::vector<mpi_cell_t> recvCOM_;
  std::vector<int> nrecvCOM_;

public:

  tree_fmm(){recvCOM_.clear(); nrecvCOM_.clear();};
  ~tree_fmm(){recvCOM_.clear(); nrecvCOM_.clear();};

  /**
   * @brief Exchange the cells required for the FMM computation
   * - Get the cells up to a determined mass 
   * - Send them to all the other processes
   * 
   * @param tree The tree topology we are working on  
   * @param maxMass The maximum mass to select the cells in the tree
   */
  void
  mpi_exchange_cells(
    tree_topology_t& tree,
    double maxMass)
  {
    int rank,size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    // Find in the tree each of the COM that will be concerned by
    // the FMM method, based on the mass first 
    std::vector<mpi_cell_t> vcells;
    std::vector<branch_t*> vbranches;

    tree.find_sub_cells_mass(
      tree.root(),
      maxMass,
      vbranches);

    vcells.resize(vbranches.size());

    int visited_entities = 0;

    #pragma omp parallel for reduction(+:visited_entities)
    for(int i = 0; i < vbranches.size(); ++i){
      assert(vbranches[i]->sub_entities()>0);
      vcells[i] = mpi_cell_t(
        vbranches[i]->get_coordinates(),
        vbranches[i]->bmin(),
        vbranches[i]->bmax(), 
        vbranches[i]->id()
      );
      visited_entities += vbranches[i]->sub_entities();
    }

    assert(visited_entities == tree.root()->sub_entities());

    nrecvCOM_.clear(); 
    // Gather the number of cells from everyone 
    nrecvCOM_.resize(size);
    std::vector<int> noffsets(size);
    nrecvCOM_[rank] = vcells.size();

    MPI_Allgather(MPI_IN_PLACE,1,MPI_INT,&nrecvCOM_[0],1,
      MPI_INT,MPI_COMM_WORLD);
  
    int totalrecv = 0;
    noffsets[0] = 0;
    for(int i=0;i<size;++i){
      totalrecv += nrecvCOM_[i];
      nrecvCOM_[i]*=sizeof(mpi_cell_t);
      if(i<size-1){
        noffsets[i+1] += nrecvCOM_[i]+noffsets[i];
      } // if
    } // for
  
    recvCOM_.clear();
    recvCOM_.resize(totalrecv);
    MPI_Allgatherv(&vcells[0],nrecvCOM_[rank],MPI_BYTE,
      &recvCOM_[0],&nrecvCOM_[0],&noffsets[0],MPI_BYTE,MPI_COMM_WORLD);
    
    if(rank == 0)
        std::cout<<"GRAV COM = ";
    for(auto v: nrecvCOM_)
      if(rank == 0)
        std::cout<<v/sizeof(mpi_cell_t)<<";";
    if(rank == 0)
        std::cout<<std::endl;
      
    // Check if mine are in the right order 
    for(size_t i=0;i<vcells.size();++i){
      assert(vcells[i].position == 
        recvCOM_[i+noffsets[rank]/sizeof(mpi_cell_t)].position);
    }// for
    MPI_Barrier(MPI_COMM_WORLD);
  } // mpi_exchange_cells

  /**
   * @brief Start the computation on the cells received from everyone
   * - Call the function tree_traversal_c2c on the cells 
   * - Use a Sink cell to store data 
   * @param tree The tree topology we are working on 
   * @param macangle The macangle for Multipole Method
   */
  void
  mpi_compute_fmm(
    tree_topology_t& tree,
    double macangle)
  {
    int size = MPI_Comm_size(MPI_COMM_WORLD,&size);
    int rank = MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    #pragma omp parallel for 
    for(auto cell=recvCOM_.begin(); cell < recvCOM_.end() ; ++cell){
      branch_t sink;
      sink.setPosition(cell->position);
      sink.set_bmax(cell->bmax);
      sink.set_bmin(cell->bmin);
      // Do the tree traversal, compute the cells data
      tree_traversal_c2c(tree,&sink,tree.root(),
          cell->fc,cell->dfcdr,cell->dfcdrdr,
          macangle,cell->ninterations);
    } // for 
  }

  // Gather the result from the other processes and add the forces 
  // Then apply to the particles below it 
  void 
  mpi_gather_cells(
    tree_topology_t& tree
    )
  { 
    int rank,size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
    int ncells = nrecvCOM_[rank]/sizeof(mpi_cell_t);

    std::vector<int> nrecv(size);
    std::vector<int> noffsets(size);
    std::vector<int> soffsets(size);
    noffsets[0] = 0;
    soffsets[0] = 0;
    std::fill(nrecv.begin(),nrecv.end(),nrecvCOM_[rank]);
    for(int i=1;i<size;++i){ 
      soffsets[i] = soffsets[i-1]+nrecvCOM_[i-1]; 
      noffsets[i] = noffsets[i-1]+nrecvCOM_[rank];
    } // for

    std::vector<mpi_cell_t> recvcells(ncells*size);
    MPI_Alltoallv(&recvCOM_[0],&nrecvCOM_[0],&soffsets[0],MPI_BYTE,
      &recvcells[0],&nrecv[0],&noffsets[0],MPI_BYTE,MPI_COMM_WORLD);

    assert((int)recvcells.size()==ncells*size);

    // Reduce the sum on the COM 
    // They are in the same order from all the others 
    for(int i=1;i<size;++i){ 
      std::cout<<"SUMMING FROM "<<i<<std::endl;
      for(int j=0;j<ncells;++j) {
        assert(recvcells[j].position ==  recvcells[i*ncells+j].position );
        assert(recvcells[j].id == recvcells[i*ncells+j].id);
        recvcells[j].fc += recvcells[i*ncells+j].fc;
        for(int k=0;k<dimension*dimension;++k){
          recvcells[j].dfcdr[k] += recvcells[i*ncells+j].dfcdr[k];
          assert(recvcells[j].dfcdr[k] < 1000);
        }
        for(int k=0;k<dimension*dimension*dimension;++k){ 
          recvcells[j].dfcdrdr[k] += recvcells[i*ncells+j].dfcdrdr[k];
          assert(recvcells[j].dfcdrdr[k] < 1000);
        } // for
        recvcells[j].ninterations += recvcells[i*ncells+j].ninterations;
      } // for 
    } // for

    // Remove to zero 
    /*for(int j=0;j<ncells;++j){ 
      for(int i = 0; i < dimension; ++i){
        if(fabs(recvcells[j].fc[i]) < 1.0e-10){
          recvcells[j].fc[i] = 0.;
        }
      }
      for(int i = 0; i < dimension*dimension; ++i){
        if(fabs(recvcells[j].dfcdr[i]) < 1.0e-10){
          recvcells[j].fc[i] = 0.;
        }
      }
      for(int i = 0; i < dimension*dimension*dimension; ++i){
        if(fabs(recvcells[j].dfcdrdr[i]) < 1.0e-10){
          recvcells[j].fc[i] = 0.;
        }
      }
    }*/

    int nbody = 0;
    // Propagate in the particles from sink 
    for(int i=0;i<ncells;++i){ 
      std::vector<body*> subparts;
      // Find the branch in the local tree with the id
      branch_t * sink =  tree.get(recvcells[i].id);
      std::cout<<"Cell =  "<<i<< " sub entities = "<< sink->sub_entities() <<" fc = "<<recvcells[i].fc<<std::endl;
      
      assert(sink!=nullptr);
      point_t pos = sink->getPosition();

      sink_traversal_c2p(tree,sink,pos,
          recvcells[i].fc,recvcells[i].dfcdr,recvcells[i].dfcdrdr,
          subparts,nbody);

      // if the bodies are non local in this area 
      assert(subparts.size()!=0);

      // I interacted with all the other particles 
      if(size == 1){
        assert(recvcells[i].ninterations+subparts.size() == tree.root()->sub_entities());
      }

      // Also do the local contribution
      for(auto bi: subparts){  
        point_t grav;
        for(auto nb: subparts){  
          double dist = flecsi::distance(bi->getPosition(),nb->getPosition());
          if(dist>0.0){  
            grav += - nb->getMass()/(dist*dist*dist)*
              (bi->getPosition()-nb->getPosition());
          } // if
        } // for
        // add in the acceleration
        bi->setAcceleration(bi->getAcceleration()+grav);
      } // for
    } // for
  } // mpi_gather_cells
  
  // Compute the acceleration due to a source branch to the sink branch 
  void 
  computeAcceleration(
    point_t sinkPosition,
    point_t sourcePosition,
    double sourceMass,
    point_t& fc, 
    //point_t& dfc,
    double* jacobi,
    double* hessian
    )
  {
    //std::cout<<"computeAcceleration"<<std::endl;
    double dist = flecsi::distance(sinkPosition,sourcePosition);
    assert(dist > 0.0);
    point_t diffPos =  sinkPosition - sourcePosition;
    //if(fabs(- sourceMass/(dist*dist*dist)*(diffPos)) > 1.0e-10)
    fc +=  -sourceMass/(dist*dist*dist)*(diffPos);
    //for(int i = 0; i < dimension; ++i){
    //  if(fabs((-sourceMass/(dist*dist*dist)*(diffPos))[i]) > 1.0e-10){
    //    fc[i] += (-sourceMass/(dist*dist*dist)*(diffPos))[i];
    //  }
    //}

    double jacobicoeff = -sourceMass/(dist*dist*dist);
    for(int i=0;i<dimension;++i){ 
      for(int j=0;j<dimension;++j){
        double valjacobian = 0.;
        if(i==j){ 
          valjacobian = jacobicoeff*(1-3*diffPos[i]*diffPos[j]/(dist*dist)); 
        }else{ 
          valjacobian = jacobicoeff*(-3*diffPos[i]*diffPos[j]/(dist*dist));
        }
        assert(!std::isnan(valjacobian));
        jacobi[i*dimension+j] += valjacobian;
      }
    }
    // Compute the Hessian matrix 
    double hessiancoeff = -3.0*sourceMass/(dist*dist*dist*dist*dist);
    for(int i=0;i<dimension;++i){
      int matrixPos = i*dimension*dimension;
      for(int j=0;j<dimension;++j){
        for(int k=0;k<dimension;++k){
          int position = matrixPos+j*dimension+k;
          //hessian[position] = 0.0;
          double firstterm = 0.0;
          if(i==j){
            firstterm += diffPos[k];
          } // if
          if(j==k){
            firstterm += diffPos[i];
          } // if
          if(k==i){
            firstterm += diffPos[j];
          } // if
          double valhessian = hessiancoeff * 
            ( 5.0/(dist*dist)*diffPos[i]*diffPos[j]*diffPos[k] - firstterm) ;

          hessian[position] += valhessian;
        } // for
      } // for
    } // for
  } // computeAcceleration

  void  
  tree_traversal_c2c(
    tree_topology_t& tree, 
    branch_t * sink, 
    branch_t * source,
    std::vector<branch_id> neighbors_branches, 
    point_t& fc, 
    double* jacobi,
    double* hessian,
    double& macangle,
    int& ninter)
  {
    std::stack<branch_t*> stk;
    stk.push(source);

    while(!stk.empty()){
      branch_t * cur = stk.top();
      stk.pop();
     
      // If the same box, stop
      // Or if inside the sink, stop 
      if((sink->bmin()==cur->bmin()&&sink->bmax()==cur->bmax()) ||
        (sink->bmin()<cur->bmin()&&sink->bmax()>cur->bmax())){
        continue;
      } // if
  
      if(geometry_t::box_MAC(
        cur->getPosition(),
        sink->getPosition(),
        cur->bmin(),
        cur->bmax(),
        macangle)){
        computeAcceleration(sink->getPosition(),cur->getPosition(),
          cur->getMass(),fc,jacobi,hessian);
        ninter+=cur->sub_entities();
        std::cout<<"Using box"<<std::endl;
      }else{
        if(cur->is_leaf()){
          for(auto bi: *cur){
            if(bi->is_local()){
              // Add particle to vector based on ID 
              
              // Check if particle is not inside my radius 
              //if(!(bi->getPosition() < sink->bmax() &&
              //  bi->getPosition() > sink->bmin())){
              //  ninter++;
              //  computeAcceleration(sink->getPosition(),bi->getPosition(),
              //    bi->getMass(),fc,jacobi,hessian);
              } // if
            } // if 
          } // for
        }else{
          for(int i=0;i<(1<<dimension);++i){
            if(tree.child(cur,i)->sub_entities() > 0){
              stk.push(tree.child(cur,i));
            } // if
          } // for
        } // if
      //} // if
    } // while
  } // tree_traversal_c2c

  void 
  sink_traversal_c2p(
    tree_topology_t& tree,
    branch_t *b,
    point_t& sinkPosition,
    point_t& fc, 
    double* jacobi,
    double* hessian,
    std::vector<body*>& neighbors,
    int& nbody
    )
  {
    std::stack<branch_t *> stk;
    stk.push(b);

    while(!stk.empty()){
      branch_t * cur = stk.top();
      stk.pop();
      if(cur->is_leaf()){
        // Apply the expansion on the bodies
        for(auto bi: *cur){
          //Do not expand on non local particles
          if(!bi->is_local()){
            continue;
          } // if
          point_t diffPos = bi->getPosition() - sinkPosition;
          point_t grav = fc;
          // The Jacobi 
          for(int i=0;i<dimension;++i){
            for(int  j=0;j<dimension;++j){
              grav[i] += jacobi[i*dimension+j]*diffPos[j];
            } // for
          } // for
          // The hessian 
          double tmpMatrix[dimension*dimension] = {};
          for(int i=0;i<dimension;++i){
            for(int j=0;j<dimension;++j){
              for(int k=0;k<dimension;++k){
                tmpMatrix[i*dimension+j] += 
                  diffPos[k]*hessian[i*dimension*dimension+j*dimension+k];
              } // for
            } // for
          } // for
          double tmpVector[dimension] = {};
          for(int i=0;i<dimension;++i){
            for(int j=0;j<dimension;++j){
              tmpVector[j] += tmpMatrix[i*dimension+j]*diffPos[i];
            } // for
          } // for
          for(int i=0;i<dimension;++i){
            grav[i] += 0.5*tmpVector[i];
          } // for
          neighbors.push_back(bi->getBody());
          bi->getBody()->setAcceleration(grav+bi->getBody()->getAcceleration());
          nbody++;
        } // for
      }else{
        for(int i=0;i<(1<<dimension);++i){
          if(tree.child(cur,i)->sub_entities() > 0){
            stk.push(tree.child(cur,i));
          } // if
        } // for
      } // else
    } // while
  } // sink_traversal_c2p
}; // class tree_fmm 

#endif // _mpisph_tree_fmm_h_
