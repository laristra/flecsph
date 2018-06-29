/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@ @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@///// /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@      /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@ /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////  /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@      /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@      /@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
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
#include "utils.h"

using namespace mpi_utils;


struct body_holder_fmm_t{
  static const size_t dimension = gdimension;
  using element_t = type_t; 
  using point_t = flecsi::point<element_t, dimension>;

  point_t position; 
  int owner; 
  double mass;
  long int id; 
  branch_id_t id_sink;
};

template<
  typename T,
  size_t D
>
class tree_fmm
{

  using geometry_t = flecsi::topology::tree_geometry<T,D>;
  static const size_t dimension = D;

  /**
   * @brief      Structure for COM communication during the gravitation
   *             computation process
  */
  struct mpi_cell_t{
    point_t position;
    double mass; 
    point_t fc = {};
    double dfcdr[dimension*dimension] = {};
    double dfcdrdr[dimension*dimension*dimension] = {};
    point_t bmax;
    point_t bmin;
    branch_id_t id;
    int ninterations = 0;
    int owner;
  
    mpi_cell_t(
      point_t position_,
      double mass_,
      point_t bmin_,
      point_t bmax_,
      branch_id_t id_,
      int owner_)
    {
      position = position_;
      mass = mass_;
      id = id_;
      bmax = bmax_;
      bmin = bmin_;
      owner = owner_;
    };

    mpi_cell_t(
      const mpi_cell_t& m1)
    {
      // Copy everythings 
      position = m1.position; 
      mass = m1.mass; 
      fc = m1.fc; 
      memcpy(dfcdr,m1.dfcdr,sizeof(double)*dimension*dimension);
      memcpy(dfcdrdr,m1.dfcdrdr,sizeof(double)*dimension*dimension*dimension);
      bmax = m1.bmax; 
      bmin = m1.bmin; 
      id = m1.id; 
      ninterations = m1.ninterations; 
      owner = m1.owner;  
    }

    mpi_cell_t& 
    operator=(
      const mpi_cell_t& m1)
    {
      // Copy everythings 
      position = m1.position; 
      mass = m1.mass; 
      fc = m1.fc; 

      memcpy(dfcdr,m1.dfcdr,sizeof(double)*dimension*dimension);
      memcpy(dfcdrdr,m1.dfcdrdr,sizeof(double)*dimension*dimension*dimension);
      
      bmax = m1.bmax; 
      bmin = m1.bmin; 
      id = m1.id; 
      ninterations = m1.ninterations; 
      owner = m1.owner;  
      return *this;
    }

    // Default 
    mpi_cell_t() {};

  };

private:

  // Used in the COM computation and exchange
  std::vector<mpi_cell_t> recvCOM_;
  std::vector<int> nrecvCOM_;

  std::vector<body_holder_fmm_t> send_particles_;
  std::vector<int> send_particles_count_;
  std::vector<body_holder_fmm_t> recv_particles_;


public:

  /**
   * @brief      Constructs the object.
   */
  tree_fmm(){recvCOM_.clear(); nrecvCOM_.clear();};

  /**
   * @brief      Destroys the object.
   */
  ~tree_fmm(){recvCOM_.clear(); nrecvCOM_.clear();};


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
    double macangle,
    double epsilon)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    double time = omp_get_wtime();

    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Reset the particles to.receive and send
    send_particles_.clear();
    recv_particles_.clear();
    send_particles_count_.clear();
    send_particles_count_.resize(size);
    std::fill(send_particles_count_.begin(),
      send_particles_count_.end(),0);

    #pragma omp parallel
    {
      std::vector<body_holder_fmm_t> local;
      std::vector<int> local_count(size);

      #pragma omp for 
      for(auto cell=recvCOM_.begin(); cell < recvCOM_.end() ; ++cell){
        branch_t sink;
        sink.setPosition(cell->position);
        sink.set_bmax(cell->bmax);
        sink.set_bmin(cell->bmin);
        sink.setMass(cell->mass);
        // Do the tree traversal, compute the cells data
        tree_traversal_c2c(tree,&sink,tree.root(),
            cell->fc,cell->dfcdr,cell->dfcdrdr,
            macangle,cell->ninterations,
            cell->owner,cell->id,epsilon,
            local_count, local);
      } // for 
     
      #pragma omp critical
      {
        // Copy in 
        send_particles_.insert(send_particles_.end(),
        local.begin(),local.end());
        // Accumulate
        for(int i = 0 ; i < size; ++i){
          send_particles_count_[i] += local_count[i];
        }
        local.clear();
        local_count.clear();
      }
    }
    
    sort(send_particles_.begin(),send_particles_.end(),
      [&rank](const body_holder_fmm_t& v1, const body_holder_fmm_t&v2){
        //assert(rank != v1.owner && rank != v2.owner);
        if(v1.owner == v2.owner){
          return v1.id < v2.id;
        }
        return v1.owner < v2.owner;
      });

    assert(send_particles_count_[rank] == 0);
    MPI_Barrier(MPI_COMM_WORLD);
    clog_one(trace)<<"mpi_compute_fmm in "<<omp_get_wtime()-time<<
      "s"<<std::endl<<std::flush;
  }

  /**
   * @brief      Calculates the gravitation between two particles. The result is
   *             store in the body bi based on position and mass of the other
   *             particle.
   *
   * @param      bi       The body sink, result is store in it 
   * @param[in]  pos_nb   The position of the source particle
   * @param[in]  mass_nb  The mass of the source particle
   */
  void computeAcceleration_direct(
    body* bi,
    const point_t& pos_nb,
    const double& mass_nb)
  {
    double dist = flecsi::distance(bi->getPosition(),pos_nb);
    if(dist > 0){
      bi->setAcceleration(bi->getAcceleration()-
        mass_nb/(dist*dist*dist)*(bi->getPosition()-pos_nb));
    }
  }
  /**
   * @brief      Compute the acceleration due to a source branch to the sink
   *             branch
   *
   * @param[in]  sinkPosition    The sink position
   * @param[in]  sourcePosition  The source position
   * @param[in]  sourceMass      The source mass
   * @param      fc              The gravitation force (function fc)
   * @param      jacobi          The jacobi matrix 
   * @param      hessian         The hessian matrix
   */
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
    double dist = flecsi::distance(sinkPosition,sourcePosition);
    assert(dist > 0.0);
    point_t diffPos =  sinkPosition - sourcePosition;
    fc +=  -sourceMass/(dist*dist*dist)*(diffPos);

    double jacobicoeff = -sourceMass/(dist*dist*dist);
    for(size_t i=0;i<dimension;++i){ 
      for(size_t j=0;j<dimension;++j){
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
    for(size_t i=0;i<dimension;++i){
      int matrixPos = i*dimension*dimension;
      for(size_t j=0;j<dimension;++j){
        for(size_t k=0;k<dimension;++k){
          int position = matrixPos+j*dimension+k;
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

  /**
   * @brief      Exchange the cells required for the FMM computation
   * - Get the cells up to a determined mass
   * - Send them to all the other processes
   *
   * @param      tree     The tree topology we are working on
   * @param      maxMass  The maximum mass to select the cells in the tree
   */
  void
  mpi_exchange_cells(
    tree_topology_t& tree,
    double maxMass)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    double time = omp_get_wtime();


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

    uint64_t visited_entities = 0;
    
    #pragma omp parallel for reduction(+:visited_entities)
    for(size_t i = 0; i < vbranches.size(); ++i){
      assert(vbranches[i]->sub_entities()>0);
      vcells[i] = mpi_cell_t(
        vbranches[i]->get_coordinates(),
        vbranches[i]->getMass(),
        vbranches[i]->bmin(),
        vbranches[i]->bmax(), 
        vbranches[i]->id(),
        rank
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
    
    clog_one(trace)<<"FMM COM = ";
    for(auto v: nrecvCOM_)
     clog_one(trace)<<v/sizeof(mpi_cell_t)<<";";
    clog_one(trace)<<std::endl;

    // Check if mine are in the right order 
    #pragma omp parallel for 
    for(size_t i=0;i<vcells.size();++i){
      assert(vcells[i].position == 
        recvCOM_[i+noffsets[rank]/sizeof(mpi_cell_t)].position);
    }// for

    MPI_Barrier(MPI_COMM_WORLD);
    clog_one(trace)<<"mpi_exchange_cells in "<<omp_get_wtime()-time<<
      "s"<<std::endl<<std::flush;
  } // mpi_exchange_cells


  /**
   * @brief      Gather the COM from other processes and the particles needed
   *             for the local computation.
   *
   * @param      tree          The tree
   * @param      macangle      The macangle
   * @param      totalnbodies  The totalnbodies
   */
  void 
  mpi_gather_cells(
    tree_topology_t& tree,
    double& macangle,
    int64_t & totalnbodies)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    double time = omp_get_wtime();

    int rank,size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Gather the distant particles 
    mpi_alltoallv(send_particles_count_,send_particles_,recv_particles_);

    // Gather the center of mass
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
      #pragma omp parallel for 
      for(int j=0;j<ncells;++j) {
        assert(recvcells[j].position ==  recvcells[i*ncells+j].position );
        assert(recvcells[j].id == recvcells[i*ncells+j].id);
        recvcells[j].fc += recvcells[i*ncells+j].fc;
        for(long unsigned int k=0;k<dimension*dimension;++k){
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

    // Propagate in the particles from sink 
    #pragma omp parallel for
    for(int i=0;i<ncells;++i){
      std::vector<body*> subparts;
      // Find the branch in the local tree with the id
      branch_t * sink =  tree.get(recvcells[i].id);
      
      assert(sink!=nullptr);
      point_t pos = sink->getPosition();

      sink_traversal_c2p(tree,sink,pos,
          recvcells[i].fc,recvcells[i].dfcdr,recvcells[i].dfcdrdr,subparts);

      // if the bodies are non local in this area 
      assert(subparts.size()!=0);
      assert(recvcells[i].ninterations <= totalnbodies);

      // Also gather the particles and proceed to the direct computation
      tree_traversal_p2p(tree,sink,tree.root(),subparts,
          recvcells[i].ninterations,macangle);

      int distant_contrib = 0;
      // Add the contribution of received particles 
      //#pragma omp parallel for reduction(+:distant_contrib)
      for(long unsigned int j = 0 ; j < recv_particles_.size(); ++j){
        if(recvcells[i].id == recv_particles_[j].id_sink){
          tree_traversal_p2p_distant(
            tree,sink,
            recv_particles_[j].position,
            recv_particles_[j].mass); 
          distant_contrib++;
        }
      }
      recvcells[i].ninterations += distant_contrib;
      if( totalnbodies != recvcells[i].ninterations ){
        std::cout<<rank<<" i="<<i<<" ntinter = "<<recvcells[i].ninterations
        <<"/"<<totalnbodies<<std::endl<<std::flush;
      }
      assert(recvcells[i].ninterations == totalnbodies);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    clog_one(trace)<<"mpi_gather_cells in "<<omp_get_wtime()-time<<
      "s"<<std::endl<<std::flush;
  } // mpi_gather_cells

private:

  /**
   * @brief      Compute interaction between the subparticles of a branch. 
   * We use the MAC angle to avoid COM already used in the c2c function.
   *
   * @param      tree      The tree
   * @param      sink      The sink
   * @param      source    The source
   * @param      subparts  The subparts, neighbors of the current particle
   * @param      ninter    The ninter
   * @param      macangle  The macangle
   */
  void 
  tree_traversal_p2p(
    tree_topology_t& tree,
    branch_t * sink,
    branch_t * source,
    std::vector<body*>& subparts,
    int& ninter,
    double& macangle)
  {
    // Tree traversal and interact only if close enough 
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    std::stack<branch_t*> stk;
    stk.push(source);

    while(!stk.empty()){
      branch_t * cur = stk.top();
      stk.pop();
  
      if(geometry_t::box_MAC(
        cur->getPosition(),
        sink->getPosition(),
        cur->bmin(),
        cur->bmax(),
        macangle)){
        // Already done 
      }else{
        if(cur->is_leaf()){
          for(auto bi: *cur){
            if(bi->is_local()){
              ninter++;
              // Apply to all subbodies of this sink
              for(auto b: subparts){
                computeAcceleration_direct(
                  b,
                  bi->getPosition(),
                  bi->getMass());
              } // for
            } // if
          } // for
        }else{
          for(int i=0;i<(1<<dimension);++i){
            if(tree.child(cur,i)->sub_entities() > 0){
              stk.push(tree.child(cur,i));
            } // if
          } // for
        } // if
      } // if
    } // while
  }

  void 
  tree_traversal_p2p_distant(
    tree_topology_t& tree,
    branch_t * sink,
    point_t& p_source,
    double& m_source)
  {
    // Tree traversal and interact only if close enough 
    std::stack<branch_t*> stk;
    stk.push(sink);

    while(!stk.empty()){
      branch_t * cur = stk.top();
      stk.pop();
  
      if(cur->is_leaf()){
        for(auto bi: *cur){
          if(bi->is_local()){
            computeAcceleration_direct(
              bi->getBody(),
              p_source,
              m_source);
          } // if
        } // for
      }else{
        for(int i=0;i<(1<<dimension);++i){
          if(tree.child(cur,i)->sub_entities() > 0){
            stk.push(tree.child(cur,i));
          } // if
        } // for
      } // if
    } // while
  }

  bool
  boundaries_intersect(
    point_t bmax1, 
    point_t bmin1, 
    point_t bmax2,
    point_t bmin2)
  { 
    point_t w1 = bmax1 - bmin1; 
    point_t w2 = bmax2 - bmin2; 

    for(int i = 0 ; i < gdimension; ++i){
      if((bmin2[i] >= bmin1[i] + w1[i]) ||
      (bmin2[i] + w2[i] <= bmin1[i])){
        return false;
      }
    }
    return true; 
  }

  void  
  tree_traversal_c2c(
    tree_topology_t& tree, 
    branch_t * sink, 
    branch_t * source,
    point_t& fc, 
    double* jacobi,
    double* hessian,
    double& macangle,
    int& ninter,
    int owner,
    branch_id_t& sink_id,
    double epsilon,
    std::vector<int>& particles_count,
    std::vector<body_holder_fmm_t>& particles)
  {
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    std::stack<branch_t*> stk;
    stk.push(source);
    // Distance between the two points
    double distance = flecsi::distance(sink->bmin(),sink->bmax())/2.;
    point_t center_sink = (sink->bmax()+sink->bmin())/2.;

    while(!stk.empty()){
      branch_t * cur = stk.top();
      stk.pop();


      if(geometry_t::box_MAC(
        cur->getPosition(),
        sink->getPosition(),
        cur->bmin(),
        cur->bmax(),
        macangle))
      {

        //point_t center_cur = (cur->bmax()+cur->bmin())/2.;

        // In every case if too close, just perform the normal computation
        //if(flecsi::distance(center_cur,center_sink) <
          //distance+epsilon && rank != owner){
          // Add the sub particles 
          //std::vector<body_holder*> sub_entities; 
          //tree.get_sub_entities_local(cur,sub_entities);
          //for(auto bi: sub_entities){
          //  particles_count[owner]++;
          //  particles.push_back(
          //    body_holder_fmm_t{
          //      bi->getPosition(),owner,bi->getMass(),bi->getId(),sink_id
          //    }
          //  ); 
         // }
        //}else{
          computeAcceleration(sink->getPosition(),cur->getPosition(),
            cur->getMass(),fc,jacobi,hessian);
          ninter+=cur->sub_entities();
        //}
      }else{
        if(cur->is_leaf()){
          for(auto bi: *cur){
            if(bi->is_local()){
              // Only send this body if not for me 
              if(owner != rank){
                particles_count[owner]++;
                particles.push_back(
                  body_holder_fmm_t{
                    bi->getPosition(),owner,bi->getMass(),bi->getId(),sink_id
                  }
                );
              }
            } // if
          } // for
        }else{
          for(int i=0;i<(1<<dimension);++i){
            if(tree.child(cur,i)->sub_entities() > 0){
              stk.push(tree.child(cur,i));
            } // if
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
    double* jacobi,
    double* hessian,
    std::vector<body*>& neighbors)
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
          for(size_t i=0;i<dimension;++i){
            for(size_t  j=0;j<dimension;++j){
              grav[i] += jacobi[i*dimension+j]*diffPos[j];
            } // for
          } // for
          // The hessian 
          double tmpMatrix[dimension*dimension] = {};
          for(size_t i=0;i<dimension;++i){
            for(size_t j=0;j<dimension;++j){
              for(size_t k=0;k<dimension;++k){
                tmpMatrix[i*dimension+j] += 
                  diffPos[k]*hessian[i*dimension*dimension+j*dimension+k];
              } // for
            } // for
          } // for
          double tmpVector[dimension] = {};
          for(size_t i=0;i<dimension;++i){
            for(size_t j=0;j<dimension;++j){
              tmpVector[j] += tmpMatrix[i*dimension+j]*diffPos[i];
            } // for
          } // for
          for(size_t i=0;i<dimension;++i){
            grav[i] += 0.5*tmpVector[i];
          } // for
          neighbors.push_back(bi->getBody());
          bi->getBody()->setAcceleration(grav+bi->getBody()->getAcceleration());
        } // for
      }else{
        for(size_t i=0;i<(1<<dimension);++i){
          if(tree.child(cur,i)->sub_entities() > 0){
            stk.push(tree.child(cur,i));
          } // if
        } // for
      } // else
    } // while
  } // sink_traversal_c2p
}; // class tree_fmm 

#endif // _mpisph_tree_fmm_h_

