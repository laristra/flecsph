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

#define COMPUTE_NEIGHBORS 1

std::ostream&
operator<<(
  std::ostream& ostr,
  const entity_key_t& id
)
{
  id.output_(ostr);
  return ostr;
}

inline 
bool
operator==(
    const point_t& p1, 
    const point_t& p2)
{
  for(size_t i=0;i<gdimension;++i)
    if(p1[i]!=p2[i])
      return false;
  return true;
}

inline 
bool
operator!=(
    const point_t& p1, 
    const point_t& p2)
{
  for(size_t i=0;i<gdimension;++i)
    if(p1[i]!=p2[i])
      return true;
  return false;
}

inline 
point_t
operator+(
    const point_t& p, 
    const double& val)
{
  point_t pr = p;
  for(size_t i=0;i<gdimension;++i)
    pr[i]+=val;
  return pr;
}

inline 
point_t
operator-(
    const point_t& p, 
    const double& val)
{
  point_t pr = p;
  for(size_t i=0;i<gdimension;++i)
    pr[i]-=val;
  return pr;
}

inline 
bool
operator<(
    const point_t& p, 
    const point_t& q)
{ 
  for(size_t i=0;i<gdimension;++i)
    if(p[i]>q[i])
      return false;
  return true;
}

inline 
bool
operator>(
    const point_t& p, 
    const point_t& q)
{  
  for(size_t i=0;i<gdimension;++i)
    if(p[i]<q[i])
      return false;
  return true;
}

inline 
point_t 
operator*(
    const point_t& p,
    const point_t& q)
{
  point_t r = p;
  for(size_t i=0;i<gdimension;++i)
    r[i] *= q[i];
  return r;
}

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

  // Used in the COM computation and exchange
  std::vector<mpi_cell_t> recvCOM;
  std::vector<int> nrecvCOM;

  const size_t noct = 256*1024; // Number of octets used for quicksort    

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
  using point_t = flecsi::point<T,dimension>;


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
 * Functions for COM and gravitation computation
 *~---------------------------------------------------------------------------*/

  // Seek for the cells that are in the mass limit
  // Send them to all the other processes 
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
  
    // Perform a tree traversal to gather the cells  
    std::function<void(branch_t*)>traverse;
    traverse = [&tree,&traverse,&vcells,&maxMass](branch_t * b){
      // Do not consider non local branches, mass is 0
      if(b->getMass() == 0.0){
        return;
      } // if
      // If this branch is a leaf or the mass is under the acceptance mass 
      if(b->is_leaf() || b->getMass() < maxMass){
        vcells.push_back(mpi_cell_t(b->getPosition(),
          b->getBMin(),b->getBMax(),b->id()));
      }else{
        for(int i=0;i<(1<<dimension);++i){
          traverse(tree.child(b,i));
        } // for
      } // if
    }; // traverse
    traverse(tree.root());
 
    nrecvCOM.clear(); 
    // Gather the number of cells from everyone 
    nrecvCOM.resize(size);
    std::vector<int> noffsets(size);
    nrecvCOM[rank] = vcells.size();
    MPI_Allgather(&nrecvCOM[rank],1,MPI_INT,&nrecvCOM[0],1,
      MPI_INT,MPI_COMM_WORLD);
  
    int totalrecv = 0;
    noffsets[0] = 0;
    for(int i=0;i<size;++i){
      totalrecv += nrecvCOM[i];
      nrecvCOM[i]*=sizeof(mpi_cell_t);
      if(i<size-1){
        noffsets[i+1] += nrecvCOM[i]+noffsets[i];
      } // if
    } // for
  
    recvCOM.clear();
    recvCOM.resize(totalrecv);
    MPI_Allgatherv(&vcells[0],nrecvCOM[rank],MPI_BYTE,
      &recvCOM[0],&nrecvCOM[0],&noffsets[0],MPI_BYTE,MPI_COMM_WORLD);
    
    // Check if mine are in the right order 
    for(size_t i=0;i<vcells.size();++i){
      assert(vcells[i].position == 
        recvCOM[i+noffsets[rank]/sizeof(mpi_cell_t)].position);
    }// for
  } // mpi_exchange_cells

  // Compute the contribution of this process on the cells sended by the 
  // other processes
  void
  mpi_compute_fmm(
    tree_topology_t& tree,
    double macangle)
  {
    #pragma omp parallel for 
    for(auto cell=recvCOM.begin(); cell < recvCOM.end() ; ++cell){
      branch_t sink;
      sink.setPosition(cell->position);
      sink.setBMax(cell->bmax);
      sink.setBMin(cell->bmin);
      //sink.setId(cell->id); 
      // Do the tree traversal, compute the cells data
      tree_traversal_c2c(tree,&sink,tree.root(),
          cell->fc,cell->dfcdr,cell->dfcdrdr,
          macangle);
    } // for 
  }

#if 0
  // Compute the ghosts particles needed for the gravity computation 
  // Then send them to the process that reqiure it
  void 
  mpi_gather_ghosts_com(
    tree_topology_t& tree,
    std::array<point_t,2>& range)
  {
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    std::vector<body_holder> sendbh; 
    std::vector<int> nsendbh(size);
    std::vector<int> soffsets(size);
    std::vector<int> nrecvbh(size);
    std::vector<int> roffsets(size);

    nsendbh[rank] =0;
    int position = 0;

    // Go through all the particles but ignore the particles on this process
    for(int i=0;i<size;++i){
      int ncells = nrecvCOM[i]/sizeof(mpi_cell_t);
      if(i!=rank){
        for(int j=0;j<ncells;++j){
          //std::cout<<vcells[position].bmin<<vcells[position].bmax<<std::endl;
          auto ents = tree.find_in_box(recvCOM[position].bmin,
              recvCOM[position].bmax);
          auto vecents = ents.to_vec();
          // Just add the ones that are local 
          for(auto bi: vecents){
            if(bi->is_local()){
              nsendbh[i]++;
              sendbh.push_back(body_holder(bi->getPosition(),
                  nullptr,rank,bi->getMass()));
            } // if
          } // for
          position++;
        } // for
      }else{
        position+=ncells;
      } // if
    } // for

    std::vector<body_holder>::iterator iter = sendbh.begin();
    // Do a sort and unique per interval 
    for(int i=0;i<size;++i){
      if(i==rank){
        assert(nsendbh[i]==0);
        continue;
      } // if
      // First sort 
      std::sort(iter,iter+nsendbh[i],
        [&range](auto& left, auto &right){ 
          return entity_key_t(range,left.getPosition()) <
            entity_key_t(range,right.getPosition());
        });
      auto last = std::unique(iter,iter+nsendbh[i],
        [&range](auto& left, auto& right){
          return entity_key_t(range,left.getPosition())
            == entity_key_t(range,right.getPosition());
        });
      int eraselements = std::distance(last,iter+nsendbh[i]);
      sendbh.erase(last,iter+nsendbh[i]);
      // Change the number of elements 
      nsendbh[i] -= eraselements;
      iter = last;
    } // for

    // Gather the send in recv 
    MPI_Alltoall(&nsendbh[0],1,MPI_INT,&nrecvbh[0],1,MPI_INT,MPI_COMM_WORLD);

    int totalrecv = 0;
    // Generate the offsets 
    for(int i=0;i<size;++i){
      totalrecv+= nrecvbh[i];
      nsendbh[i]*=sizeof(body_holder);
      nrecvbh[i]*=sizeof(body_holder);
      if(i>1){ 
        soffsets[i] = soffsets[i-1]+nsendbh[i-1];
        roffsets[i] = roffsets[i-1]+nrecvbh[i-1];
      } // if
    } // for

    std::vector<body_holder> recvbh(totalrecv);
    MPI_Alltoallv(&sendbh[0],&nsendbh[0],&soffsets[0],MPI_BYTE,
      &recvbh[0],&nrecvbh[0],&roffsets[0],MPI_BYTE,MPI_COMM_WORLD);

#ifdef OUTPUT
    std::cout<<rank<<": gathered="<<recvbh.size()<<std::endl;
#endif
    // Create a local tree 
    tree_topology_t localtree(range[0],range[1]);
    // Add in the tree 
    for(auto bi: recvbh){
      auto nbi = localtree.make_entity(bi.getPosition(),nullptr,bi.getOwner(),
        bi.getMass());
      localtree.insert(nbi);
    } // for

    int ncells = nrecvCOM[rank]/sizeof(mpi_cell_t);
    for(int i=0;i<ncells;++i){
      auto subents = localtree.find_in_box(recvCOM[i].bmin,recvCOM[i].bmax);
      auto subentsapply = tree.find_in_box(recvCOM[i].bmin,recvCOM[i].bmax);
      auto subvec = subents.to_vec();
      auto subvecapply = subentsapply.to_vec();
      for(auto bi: subvecapply){  
        if(bi->is_local( )){
          point_t grav = bi->getBody()->getGravForce();
          for(auto nb: subvec){  
            double dist = flecsi::distance(bi->getPosition(),nb->getPosition());
            if(dist>0.0){  
              grav += - nb->getMass()/(dist*dist*dist)*
                (bi->getPosition()-nb->getPosition());
            } // if
          } // for
          bi->getBody()->setGravForce(grav);
        } // if
      } // for
    } // for
  } // mpi_gather_ghosts_com
#endif

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
  
    int totalrecv = nrecvCOM[rank]*size;
    int ncells = nrecvCOM[rank]/sizeof(mpi_cell_t);

    std::vector<int> nrecv(size);
    std::vector<int> noffsets(size);
    std::vector<int> soffsets(size);
    noffsets[0] = 0;
    soffsets[0] = 0;
    std::fill(nrecv.begin(),nrecv.end(),nrecvCOM[rank]);
    for(int i=1;i<size;++i){ 
      soffsets[i] = soffsets[i-1]+nrecvCOM[i-1]; 
      noffsets[i] = noffsets[i-1]+nrecvCOM[rank];
    } // for

    std::vector<mpi_cell_t> recvcells(ncells*size);
    MPI_Alltoallv(&recvCOM[0],&nrecvCOM[0],&soffsets[0],MPI_BYTE,
      &recvcells[0],&nrecv[0],&noffsets[0],MPI_BYTE,MPI_COMM_WORLD);

    assert((int)recvcells.size()==ncells*size);
    //std::cout<<rank<<": receiv total="<<recvcells.size()<<std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);

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
  
    //std::cout<<rank<<": updated"<<std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);
  
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
      assert(subparts.size()!=0);
      // Also apply the subcells
      for(auto bi: subparts){  
        point_t grav = point_t{};
        for(auto nb: subparts){  
          double dist = flecsi::distance(bi->getPosition(),nb->getPosition());
          if(dist>0.0){  
            grav += - nb->getMass()/(dist*dist*dist)*
              (bi->getPosition()-nb->getPosition());
          } // if
        } // for
        //bi->setGravForce(grav);
        // add in the acceleration
        bi->setAcceleration(bi->getAcceleration()+grav);
      } // for
    } // for
    //std::cout<<rank<<": computed"<<std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);
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

  bool 
  box_intersection(
    point_t& sinkBMin,
    point_t& sinkBMax,
    point_t& sourceBMin, 
    point_t& sourceBMax)
  {
    return (sinkBMin[0]<=sourceBMax[0]&&sinkBMax[0]>=sourceBMin[0])&&
           (sinkBMin[1]<=sourceBMax[1]&&sinkBMax[1]>=sourceBMin[1])&&
           (sinkBMin[2]<=sourceBMax[2]&&sinkBMax[2]>=sourceBMin[2]);
  }

  bool
  MAC(
    branch_t * sink,
    branch_t * source,
    double macangle)
  {
    double dmax = flecsi::distance(source->getBMin(),source->getBMax());
    double disttoc = flecsi::distance(
        sink->getPosition(),source->getPosition());
    return dmax/disttoc < macangle;
  }

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
    // Check if it is the same key 
    if(sink->id()==source->id())
    {
      std::cout<<"Same Id"<<std::endl;
    }
    if(source->getMass() == 0.0){
      return;
    } // if
   
    // If the same box, stop
    if(sink->getBMin()==source->getBMin()&&sink->getBMax()==source->getBMax()){
      return;
    } // if

    // If inside the sink, stop 
    if(sink->getBMin()<source->getBMin()&&sink->getBMax()>source->getBMax()){
      return;
    } // if
  
    if(MAC(sink,source,macangle)){
      computeAcceleration(sink->getPosition(),source->getPosition(),
        source->getMass(),fc,jacobi,hessian);
    }else{
      if(source->is_leaf()){
        for(auto bi: *source){
          if(bi->is_local()){
            // Check if particle is inside my radius 
            if(!(bi->getPosition() < sink->getBMax() &&
              bi->getPosition() > sink->getBMin())){
              computeAcceleration(sink->getPosition(),bi->getPosition(),
                bi->getMass(),fc,jacobi,hessian);
            } // if
          } // if
        } // for
      }else{
        for(int i=0;i<(1<<dimension);++i){
          tree_traversal_c2c(tree,sink,tree.child(source,i),
            fc,jacobi,hessian,macangle);
        } // for
      } // if
    } // if
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
    if(b->getMass() <= 0.0){
      return;
    } // if
    if(b->is_leaf()){
      // Apply the expansion on the bodies
      for(auto bi: *b){
        if(!bi->is_local()){
          continue;
        } // if
        //point_t diffPos = sinkPosition-bi->getPosition();
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
        //bi->getBody()->setGravForce(grav); 
        // add this contribution to acceleration
        bi->getBody()->setAcceleration(grav+bi->getBody()->getAcceleration());
        
        nbody++;
      } // for
    }else{
      for(int i=0;i<(1<<dimension);++i){
        sink_traversal_c2p(tree,tree.child(b,i),sinkPosition,fc,jacobi,hessian,
          neighbors,nbody);
      } // for
    } // if
  } // sink_traversal_c2p


  // Compute the center of mass values from the particles 
  void tree_traversal_com(
      tree_topology_t& tree)
  {
    std::function<void(branch_t*)>traverse;
    // Traversal function using DFS
    traverse = [&tree,&traverse](branch_t * b){
      double mass = 0.0;
      point_t com = {};
      point_t bmax = {-99999,-99999,-99999};
      point_t bmin = {99999,99999,99999};
      if(b->is_leaf()){
        for(auto child: *b){
          // Only for local particles 
          if(child->is_local()){
            double h = child->getBody()->getSmoothinglength();
            assert(child->getMass()>0.0);
            com += child->getMass()*child->getPosition();
            mass += child->getMass();
            for(size_t i=0;i<dimension;++i){
              if(bmax[i] < child->getPosition()[i]){
                bmax[i] = child->getPosition()[i];
              } // if
              if(bmin[i] > child->getPosition()[i]){
                bmin[i] = child->getPosition()[i];
              } // if
            } // for
          } // if
        } // for 
        if(mass > 0.0){
          com = com / mass;
        } // if 
      }else{
        for(int i=0;i<(1<<dimension);++i){
          auto branch = tree.child(b,i);
          traverse(branch);
          // Add this children position and coordinates
          com += branch->getMass()*branch->getPosition();
          mass += branch->getMass();
          for(size_t i=0;i<dimension;++i){
            if(bmax[i] < branch->getBMax()[i]){
              bmax[i] = branch->getBMax()[i];
            } // if 
            if(bmin[i] > branch->getBMin()[i]){
              bmin[i] = branch->getBMin()[i];
            } // if
          } // for
        } // for
        if(mass > 0.0){
          com = com / mass;
        } // if
      } // if
      b->setMass(mass);
      assert(!std::isnan(mass));
      assert(mass>=0.0);
      b->setPosition(com);
      for(size_t i=0;i<dimension;++i){
        assert(!std::isnan(com[i]));
      } // for
      b->setBMax(bmax);
      b->setBMin(bmin);
    }; // traverse function
    
    traverse(tree.root());
    
#ifdef DEBUG
    // For Debug, check if mass conserved 
    static double checkmass = 0.0;
    double mass = tree.root()->getMass();
    MPI_Allreduce(MPI_IN_PLACE,&mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(checkmass == 0.0){
      checkmass = mass;
    }else{
      assert(fabs(checkmass - mass) < 1.0e-15 );
    }
#endif 
  } // tree_traversal_com

}; // class tree_colorer 

#endif // _mpisph_tree_colorer_h_
