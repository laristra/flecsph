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

/**
 * @brief      Structure for COM communication during the gravitation
 *             computation process
 */
struct mpi_cell_t{
  point_t position;
  point_t fc = {};
  double dfcdr[9] = {};
  double dfcdrdr[27] = {};
  point_t bmax;
  point_t bmin;
  branch_id_t id;
  
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

template<
  typename T,
  size_t D
>
class tree_fmm
{

  using point_t = flecsi::point<T,D>;
  using geometry_t = flecsi::topology::tree_geometry<T,D>;
  static const size_t dimension = D;

private:

  // Used in the COM computation and exchange
  std::vector<mpi_cell_t> recvCOM_;
  std::vector<int> nrecvCOM_;
  double mass_criterion_ = 1.0e-40;

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
      mass_criterion_,
      vbranches);

    vcells.resize(vbranches.size());

    #pragma omp parallel for 
    for(int i = 0; i < vbranches.size(); ++i){
      assert(vbranches[i]->sub_entities()>0);
      vcells[i] = mpi_cell_t(
        vbranches[i]->get_coordinates(),
        vbranches[i]->bmin(),
        vbranches[i]->bmax(), 
        vbranches[i]->id()
      );
    }

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
    int size = MPI_Comm_rank(MPI_COMM_WORLD,&size);
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
          macangle);
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
  
    int totalrecv = nrecvCOM_[rank]*size;
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
      for(int j=0;j<ncells;++j) {
        assert(recvcells[j].position ==  recvcells[i*ncells+j].position );
        assert(recvcells[j].id == recvcells[i*ncells+j].id);
        recvcells[j].fc += recvcells[i*ncells+j].fc;
        for(int k=0;k<27;++k){ 
          if(k<9){ 
            recvcells[j].dfcdr[k] += recvcells[i*ncells+j].dfcdr[k];
          } // if
          recvcells[j].dfcdrdr[k] += recvcells[i*ncells+j].dfcdrdr[k];
          // Check if cells are not too high 
          assert(recvcells[j].dfcdrdr[k] < 1000);
        } // for
      } // for 
    } // for
  
    int nbody = 0;
    // Propagate in the particles from sink 
    for(int i=0;i<ncells;++i){ 
      std::vector<body*> subparts;
      // Find the branch in the local tree with the id
      branch_t * sink =  tree.get(recvcells[i].id);
      assert(sink!=nullptr);
      point_t pos = sink->getPosition();

      sink_traversal_c2p(tree,sink,pos,
          recvcells[i].fc,recvcells[i].dfcdr,recvcells[i].dfcdrdr,
          subparts,nbody);

      // if the bodies are non local in this area 
      assert(subparts.size()!=0);

      // Also do the local contribution
      for(auto bi: subparts){  
        point_t grav = point_t{};
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
    fc +=  - sourceMass/(dist*dist*dist)*(diffPos);
    double jacobicoeff = -sourceMass/(dist*dist*dist);
    for(int i=0;i<3;++i){ 
      for(int j=0;j<3;++j){
        if(i==j){ 
          jacobi[i*3+j] += jacobicoeff* 
            (1 - 3*(diffPos[i])*(diffPos[j])/(dist*dist)); 
        }else{ 
          jacobi[i*3+j] += jacobicoeff*
          (- 3*(diffPos[i])*(diffPos[j])/(dist*dist));
        }
        assert(!std::isnan(jacobi[i*3+j]));
      }
    }
    // Compute the Hessian matrix 
    double hessiancoeff = -3.0*sourceMass/(dist*dist*dist*dist*dist);
    for(int i=0;i<3;++i){
      int matrixPos = i*9;
      for(int j=0;j<3;++j){
        for(int k=0;k<3;++k){
          int position = matrixPos+j*3+k;
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
          if(!(i==j==k)){
            firstterm*=3.0;
          } // if
          hessian[position] += hessiancoeff * firstterm + 
            hessiancoeff * -5.0/(dist*dist)*
            (diffPos[i])*(diffPos[j])*(diffPos[k]);
        } // for
      } // for
    } // for
  } // computeAcceleration

  void  
  tree_traversal_c2c(
    tree_topology_t& tree, 
    branch_t * sink, 
    branch_t * source, 
    point_t& fc, 
    //point_t& dfc,
    double* jacobi,
    double* hessian,
    double& macangle)
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
        //std::cout<<"Same box"<<std::endl;
        return;
      } // if
  
      if(geometry_t::box_MAC(
        cur->getPosition(),
        sink->getPosition(),
        cur->bmin(),
        cur->bmax(),
        macangle)){
        computeAcceleration(sink->getPosition(),cur->getPosition(),
          cur->getMass(),fc,jacobi,hessian);
      }else{
        if(cur->is_leaf()){
          for(auto bi: *cur){
            if(bi->is_local()){
              // Check if particle is not inside my radius 
              if(!(bi->getPosition() < sink->bmax() &&
                bi->getPosition() > sink->bmin())){
                computeAcceleration(sink->getPosition(),bi->getPosition(),
                  bi->getMass(),fc,jacobi,hessian);
              } // if
            } // if
          } // for
        }else{
          for(int i=0;i<(1<<dimension);++i){
            if(tree.child(cur,i)->sub_entities() > 0){
              stk.push(tree.child(cur,i));
            }
          } // for
        } // if
      } // if
    } // while
  } // tree_traversal_c2c

  void 
  sink_traversal_c2p(
    tree_topology_t& tree,
    branch_t *b,
    point_t& sinkPosition,
    point_t& fc, 
    //point_t& dfc, 
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
          for(int i=0;i<3;++i){
            for(int  j=0;j<3;++j){
              grav[i] += (jacobi[i*3+j]*(diffPos[j]));
            } // for
          } // for
          // The hessian 
          double tmpMatrix[9] = {};
          for(int i=0;i<3;++i){
            for(int j=0;j<3;++j){
              for(int k=0;k<3;++k){
                tmpMatrix[i*3+j] += diffPos[k]*hessian[i*9+j+k*3];
              } // for
            } // for
          } // for
          double tmpVector[3] = {};
          for(int i=0;i<3;++i){
            for(int j=0;j<3;++j){
              tmpVector[i] += tmpMatrix[i*3+j]*diffPos[j];
            } // for
          } // for
          for(int i=0;i<3;++i){
            grav[i] += 1./2.*tmpVector[i];
          } // for
          neighbors.push_back(bi->getBody());
          bi->getBody()->setAcceleration(grav+bi->getBody()->getAcceleration());
          nbody++;
        } // for
      }else{
        for(int i=0;i<(1<<dimension);++i){
          if(tree.child(cur,i)->sub_entities() > 0){
            stk.push(tree.child(cur,i));
          }
        } // for
      } // else
    } // while
  } // sink_traversal_c2p
}; // class tree_fmm 

#endif // _mpisph_tree_fmm_h_
