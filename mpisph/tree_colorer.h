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

#include "tree.h"

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
  std::vector<int>  sholders;
  std::vector<int>  soffsets;
  std::vector<int>  rholders;
  std::vector<int>  roffsets;
  std::vector<body_holder*> recvholders;
  std::vector<std::set<body_holder*>> sendholders;
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

  const size_t noct = 2048*1024; // Number of octets used for quicksort    

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
    ghosts_data.sholders.resize(size);
    ghosts_data.rholders.resize(size);
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
    ghosts_data.sholders.resize(size);
    ghosts_data.rholders.resize(size);
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
    int totalnbodies
  )
  {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

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
      return;
    } // if
    
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<std::pair<entity_key_t,int64_t>> splitters;
    generate_splitters(splitters,rbodies,totalnbodies); 

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
    if(rank == 0){
      std::cout<<"Repartition: ";
      for(auto num: totalprocbodies)
        std::cout<<num<<";";
      std::cout<<std::endl;
    } // if
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
   * @param[in]  smoothinglength  The smoothinglength, the max value of everyone
   */
  void
  mpi_branches_exchange(
    tree_topology_t& tree,
    std::vector<std::pair<entity_key_t,body>>& rbodies,
    std::vector<std::array<point_t,2>>& ranges,
    std::array<point_t,2>& range_total,
    double smoothinglength)
  {
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    
#ifdef OUTPUT
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<"Branches repartition" << std::flush;
#endif

   
    // Search for my min and max posititon
    std::array<point_t,2> range;
    local_range(rbodies,range);

    // Add the smoothing length to the max and min to find the real boundaries
    // Consider that on the edge of the box, and the sphere radius. 
    double dist_sphere_box = sqrt(2)*2*smoothinglength;
    range[0] = range[0]-dist_sphere_box;
    range[1] = range[1]+dist_sphere_box;

    // Gather the keys of everyone 
    // If it is the first time, allocate the ranges 
    if(ranges.size() == 0){
      ranges.resize(size);
    }
    
    MPI_Allgather(&range,sizeof(std::array<point_t,2>),
      MPI_BYTE,&ranges[0],sizeof(std::array<point_t,2>),
      MPI_BYTE,MPI_COMM_WORLD);

    // Now generate the sendbuffer, ordered by processes
    // for the holders 
    reset_buffers();
    std::vector<body_holder_mpi_t> sendbuffer;
    scount[rank]=0;

    // Search in the tree for each processes 
    for(int i=0;i<size;++i)
    {
      if(i==rank)
        continue;
      auto ents = tree.find_in_box_b(ranges[i][0],ranges[i][1]);
      scount[i] = ents.size();
      for(auto ent: ents){
        sendbuffer.push_back(body_holder_mpi_t{
            ent->getPosition(),
            rank,
            ent->getMass(),
	          ent->getId()});
      }
    }

    std::vector<body_holder_mpi_t> recvbuffer;
    mpi_alltoallv(scount,sendbuffer,recvbuffer); 

    // Add them in the tree 
    for(auto bi: recvbuffer)
    {
      assert(bi.owner!=rank);
      assert(bi.mass!=0.);
      auto nbi = tree.make_entity(bi.position,nullptr,bi.owner,bi.mass,bi.id);
      tree.insert(nbi);
    }

#ifdef OUTPUT
    if(rank==0)
      std::cout<<".done"<<std::endl;
#endif

  }

  /**
   * @brief      Update the local ghosts data. This function is always use 
   * after one call of mpi_compute_ghosts. It can be use several time like:
   * mpi_compute_ghosts 
   * mpi_refresh_ghosts 
   * mpi_refresh_ghosts 
   * As long as the particles did not move. 
   *
   * @param      tree   The tree
   * @param      range  The range
   */
  void mpi_refresh_ghosts(
    tree_topology_t& tree/*,*/ 
    /*std::array<point_t,2>& range*/)
  {
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
  
#ifdef OUTPUT
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<"Refresh Ghosts" << std::flush;
#endif 

    // Refresh the sendbodies with new data
    auto itsb = ghosts_data.sbodies.begin();
    for(auto proc: ghosts_data.sendholders)
    {
      for(auto bi: proc)
      {
        assert(bi->getBody()!=nullptr);
        *itsb = *(bi->getBody());
        itsb++;
      }
    }  

    MPI_Alltoallv(&ghosts_data.sbodies[0],&ghosts_data.sholders[0],
      &ghosts_data.soffsets[0],MPI_BYTE,
      &ghosts_data.rbodies[0],&ghosts_data.rholders[0],
      &ghosts_data.roffsets[0],MPI_BYTE,
      MPI_COMM_WORLD);

    // Sort the bodies based on key and or position
    std::sort(ghosts_data.rbodies.begin(),
            ghosts_data.rbodies.end(),
      [](auto& left, auto& right){
        if(entity_key_t(left.coordinates())<
          entity_key_t(right.coordinates())){
          return true;
        }
        if(entity_key_t(left.coordinates())==
          entity_key_t(right.coordinates())){
          return left.getId()<right.getId();
        }
        return false;
      });
 

    assert(ghosts_data.rbodies.size() == ghosts_data.recvholders.size()); 

    // Then link the holders with these bodies
    auto it = ghosts_data.rbodies.begin();
    if(size==1)
      assert(ghosts_data.recvholders.size() == 0); 
    for(auto& bi: ghosts_data.recvholders)
    {
      auto bh = tree.get(bi->id());
      bh->setBody(&(*it));

      assert(!bh->is_local() && "Non local particle");
      assert(entity_key_t(bh->coordinates()) == 
          entity_key_t(bh->getBody()->coordinates()) && 
          "Key different than holder");
      assert(bh->coordinates()==bh->getBody()->coordinates() &&
          "Position different than holder");
      assert(bh->getId() == bh->getBody()->getId()&& 
          "Id different from holder");     
      
      // Replace the body_holder by the received
      bh->getBody()->setPosition(bh->getPosition());
      bh->setPosition(bh->getBody()->getPosition());
      ++it;
    }   
  
#ifdef OUTPUT
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done" << std::endl << std::flush;
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
   * @param[in]  smoothinglength  The smoothinglength
   * @param      range            The range
   */
  void 
  mpi_compute_ghosts(
    tree_topology_t& tree,
    double smoothinglength
  )
  {    
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
  
#ifdef OUTPUT
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<"Compute Ghosts" << std::flush;
#endif

    // Clean the structure 
    ghosts_data.sbodies.clear();
    ghosts_data.rbodies.clear();
    ghosts_data.recvholders.clear();
    ghosts_data.sendholders.clear();
        
    assert(ghosts_data.sholders.size()==(size_t)size);

    // 1. for each processes, bucket of send, bucket of receiv 
    ghosts_data.sendholders.resize(size);

    std::vector<std::set<body_holder*>> recvholders(size);
    
    // TODO add a reduction over h 
    int64_t totalrecvbodies = 0L; 
    int64_t totalsendbodies = 0L; 
    auto treeents = tree.entities().to_vec(); 
    for(auto bi: treeents)
    {  
      if(bi->is_local())
      {
        assert(bi->getOwner() == rank);
        assert(bi->getBody()->getSmoothinglength() > 0); 
        auto bodiesneighbs = tree.find_in_radius_b(
            bi->coordinates(),
            2.*bi->getBody()->getSmoothinglength());
        assert(bodiesneighbs.size() > 0);
        for(auto nb: bodiesneighbs)
	      {
          if(!nb->is_local())
          {
	          assert(nb->getOwner()!=rank && nb->getOwner() != -1); 
            ghosts_data.sendholders[nb->getOwner()].insert(bi);
            recvholders[nb->getOwner()].insert(nb);
          } // if
        } // for
      } // if 
    } // for

    for(int i=0;i<size;++i){
      ghosts_data.sholders[i] = ghosts_data.sendholders[i].size();
      assert(ghosts_data.sholders[i]>=0);
      totalsendbodies += ghosts_data.sholders[i]; 
      ghosts_data.rholders[i] = recvholders[i].size();
      assert(ghosts_data.rholders[i]>=0);
      totalrecvbodies += ghosts_data.rholders[i]; 
    } // for
  
    // Make a vector with the recvholsters to be able to connect the pointer
    // at the end of the communication
    for(auto proc: recvholders){
      ghosts_data.recvholders.insert(ghosts_data.recvholders.end(),proc.begin(),
        proc.end());
    } // for

    // Now gather the bodies data to send in a vector 
    // Take the holders in the order 0 to n processes 
    ghosts_data.sbodies.resize(totalsendbodies);

    // Prepare offsets for alltoallv
    ghosts_data.roffsets[0]=0;
    ghosts_data.soffsets[0]=0;

    for(int i=1;i<size;++i){
      ghosts_data.roffsets[i] = ghosts_data.rholders[i-1]+
        ghosts_data.roffsets[i-1];
      ghosts_data.soffsets[i] = ghosts_data.sholders[i-1]+
        ghosts_data.soffsets[i-1]; 
    }

    ghosts_data.rbodies.resize(totalrecvbodies);

    // Convert the offsets to byte
    for(int i=0;i<size;++i){
      ghosts_data.sholders[i]*=sizeof(body);
      ghosts_data.rholders[i]*=sizeof(body);
      ghosts_data.soffsets[i]*=sizeof(body);
      ghosts_data.roffsets[i]*=sizeof(body);
    }
  
    // Sort the holders once
    std::sort(ghosts_data.recvholders.begin(),
            ghosts_data.recvholders.end(),
      [](auto& left, auto& right){
        if(entity_key_t(left->coordinates())<
          entity_key_t(right->coordinates())){
          return true;
        }
        if(entity_key_t(left->coordinates())==
          entity_key_t(right->coordinates())){          
 	        return left->getId()<right->getId();
        }
        return false;
      });
      
#ifdef OUTPUT 
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<".done"<<std::endl;
#endif
  }

/*~---------------------------------------------------------------------------*
 * Utils functions
 *~---------------------------------------------------------------------------*/


/**
 * @brief      Compute the local range of particles
 * range 0 = min range 1 = max
 *
 * @param      bodies  The bodies
 * @param      range   The range
 */
  void 
  local_range(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    std::array<point_t,2>& range)
  {
    
    range[1] = bodies.back().second.coordinates();
    range[0] = bodies.back().second.coordinates();
    
    for(auto bi: bodies){
      for(size_t i=0;i<dimension;++i){
        if(bi.second.coordinates()[i]>range[1][i])
          range[1][i] = bi.second.coordinates()[i];
        if(bi.second.coordinates()[i]<range[0][i])
          range[0][i] = bi.second.coordinates()[i];
      }
    }
  }

  /**
   * @brief      Compute the global range of all the particle system 
   *
   * @param      bodies           The bodies, local of this process
   * @param      range            The range computed in that function  
   * @param[in]  smoothinglength  The smoothinglength, the biggest of the system
   */
  void 
  mpi_compute_range(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    std::array<point_t,2>& range,
    double smoothinglength)
  {
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Compute the local range 
    std::array<point_t,2> lrange; 
    local_range(bodies,lrange);
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
 
    point_t minposition; 
    point_t maxposition; 

    for(size_t i=0;i<dimension;++i){
      minposition[i] = min[i]-2*smoothinglength;
      maxposition[i] = max[i]+2*smoothinglength;
    }

#ifdef OUTPUT
    if(rank==0)
      std::cout <<"boundaries: "<< minposition << maxposition << std::endl;
#endif 

    range[0] = minposition;
    range[1] = maxposition;
  }


/**
 * @brief      Simple version of all to all 
 * Use to generate the offsets and do the pre exchange 
 * Then realise the MPI_Alltoall call 
 *
 * @param[in]  sendcount   The sendcount
 * @param      sendbuffer  The sendbuffer
 * @param      recvbuffer  The recvbuffer
 *
 * @tparam     M           The type of data sent 
 */
  template<
    typename M>
  void 
  mpi_alltoallv(
      std::vector<int> sendcount,
      std::vector<M>& sendbuffer,
      std::vector<M>& recvbuffer
    )
  {
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size); 

    std::vector<int> recvcount(size); 
    std::vector<int> recvoffsets(size); 
    std::vector<int> sendoffsets(size); 

    // Exchange the send count 
    MPI_Alltoall(&sendcount[0],1,MPI_INT,&recvcount[0],1,MPI_INT,
      MPI_COMM_WORLD);
    
    // Generate the send and recv offsets 
    std::partial_sum(recvcount.begin(),recvcount.end(),&recvoffsets[0]); 
    // As we need an exscan, add a zero
    recvoffsets.insert(recvoffsets.begin(),0);
    
    // Then send offsets
    std::partial_sum(sendcount.begin(),sendcount.end(),&sendoffsets[0]);
    // As we need an exscan, add a zero
    sendoffsets.insert(sendoffsets.begin(),0);
    
    // Set the recvbuffer to the right size
    recvbuffer.resize(recvoffsets.back()); 
    
    // Trnaform the offsets for bytes 
    for(int i=0;i<size;++i){
      sendcount[i] *= sizeof(M);
      recvcount[i] *= sizeof(M);
      sendoffsets[i] *= sizeof(M);
      recvoffsets[i] *= sizeof(M);
    } // for
    
    // Use this array for the global buckets communication
    MPI_Alltoallv(&sendbuffer[0],&sendcount[0],&sendoffsets[0],MPI_BYTE,
      &recvbuffer[0],&recvcount[0],&recvoffsets[0],MPI_BYTE,MPI_COMM_WORLD);
  }

/**
 * @brief      Use in mpi_qsort to generate the splitters to sort the particles 
 * in the quick sort algorithm 
 * In this function we take some samplers of the total particles and the root 
 * determines the splitters
 *
 * @param      splitters  The splitters used in the qsort in mpi_qsort
 * @param[in]  rbodies  The local bodies of the process
 */
  void
  generate_splitters(
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
    size_t nsample = noct / sizeof(std::pair<entity_key_t,int64_t>);
    if(rbodies.size()<nsample){nsample = rbodies.size();}
    int chuncksize = rbodies.size()/nsample;
    
    for(size_t i=0;i<nsample;++i){
      keys_sample.push_back(std::make_pair(rbodies[chuncksize*i].first,
      rbodies[chuncksize*i].second.getId()));
    } // for
    assert(keys_sample.size()==nsample);

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
      master_nkeys = noct/sizeof(std::pair<entity_key_t,int64_t>)*size;
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
      chuncksize = master_nkeys/size;
      for(int i=0;i<size-1;++i){
        splitters[i] = master_keys[(i+1)*chuncksize];
        //std::cout<<splitters[i].first<<std::endl;
      } // for
      // Last key 
      //std::cout<<entity_key_t::last_key()<<std::endl;
    } // if

    // Bradcast the splitters 
    MPI_Bcast(&splitters[0],(size-1)*sizeof(std::pair<entity_key_t,int64_t>)
    ,MPI_BYTE,0,MPI_COMM_WORLD);
  }

  /**
   * @brief      Export to a file the current tree in memory 
   * This is useful for small number of particles to help representing the tree 
   *
   * @param      tree   The tree to output
   * @param      range  The range of the particles, use to construct entity_keys
   */
  void mpi_tree_traversal_graphviz(
    tree_topology_t & tree/*,*/
    /*std::array<point_t,2>& range*/)
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    char fname[64];
    sprintf(fname,"output_graphviz_%d.gv",rank);
    std::ofstream output;
    output.open(fname);
    output<<"digraph G {"<<std::endl;

    std::stack<branch_t*> stk;
    // Get root
    auto rt = tree.root();
    stk.push(rt);

    while(!stk.empty()){
      branch_t* cur = stk.top();
      stk.pop();
      if(!cur->is_leaf()){
        // Add the child to the stack and add for display 
        for(size_t i=0;i<(1<<dimension);++i)
        {
          auto br = tree.child(cur,i);
          stk.push(br);
          if(dimension == 3)
          {
            output<<std::oct<<cur->id().value_()
              <<"->"<<br->id().value_()<<std::dec<<std::endl;
          }else if(dimension == 1)
          {
            output<<std::bitset<64>(cur->id().value_())<<"->"<<
              std::bitset<64>(br->id().value_())<<std::endl;
          }
        }
      }else{
        for(auto ent: *cur)
        {
          entity_key_t key(/*range,*/ent->coordinates());
          output<<std::bitset<64>(cur->id().value_())<<
            "->"<<key<<std::endl;
          switch (ent->getLocality())
          {
            case 2:
              output<<key<<"[shape=box,color=blue]"<<std::endl;
              break;
            case 3:
              output<<key<<" [shape=box,color=red]"<<std::endl;
              //fprintf(output,"\"%lo\" [shape=box,color=red]\n",
              //  key.truncate_value(17));
              break;
            case 1:
              output<<key<<" [shape=box,color=green]"<<std::endl;
              //fprintf(output,"\"%lo\" [shape=box,color=green]\n",
              //  key.truncate_value(17));
              break;
            default:
              output<<key<<" [shape=circle,color=black]"<<std::endl;
              //fprintf(output,"\"%lo\" [shape=circle,color=black]\n",
              //  key.truncate_value(17));
              break;
          }
        }
      } 
    }
    output<<"}"<<std::endl;
    output.close();
  }


}; // class tree_colorer 

#endif // _mpisph_tree_colorer_h_
