/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 
#ifndef _mpisph_body_system_h_
#define _mpisph_body_system_h_

#include "tree_colorer.h"
#include "io.h"

template<
  typename T,
  size_t D
  >
class body_system{

using point_t = flecsi::point<T,D>;

public:
  body_system():totalnbodies_(0L),localnbodies_(0L),macangle_(0.0),
  maxmasscell_(1.0e-40),tree_(nullptr)
  {};
  ~body_system(){
    if(tree_ != nullptr){
      delete tree_;
    }
  };


  void setMaxmasscell(double maxmasscell){maxmasscell_ = maxmasscell;};
  void setMacangle(double macangle){macangle_ = macangle;};

  void read_bodies(
      const char * filename,
      int startiteration)
  {
    io::inputDataHDF5(localbodies_,filename,totalnbodies_,localnbodies_,
        startiteration);
  }

  void read_bodies_txt(
      const char* filename)
  {
    io::inputDataTxtRange(localbodies_,localnbodies_,totalnbodies_,filename);
  }

  void write_bodies(
      const char * filename, 
      int iter,
      bool do_diff_files = false)
  {
    io::outputDataHDF5(localbodies_,filename,iter,do_diff_files);
  }

  double getSmoothinglength()
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Choose the smoothing length to be the biggest from everyone 
    double smoothinglength = 0;
    for(auto bi: localbodies_){
      if(smoothinglength < bi.second.getSmoothinglength()){
        smoothinglength = bi.second.getSmoothinglength();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&smoothinglength,1,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);

    if(rank==0){
      std::cout<<"H="<<smoothinglength<<std::endl;
    }
    return smoothinglength;

  }

  std::array<point_t,2>& getRange()
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Destroy the previous tree
    if(tree_ !=  nullptr){
      delete tree_;
    }
    
  
    // Choose the smoothing length to be the biggest from everyone 
    double smoothinglength = 0;
    for(auto bi: localbodies_){
      if(smoothinglength < bi.second.getSmoothinglength()){
        smoothinglength = bi.second.getSmoothinglength();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&smoothinglength,1,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);

    if(rank==0){
      std::cout<<"H="<<smoothinglength<<std::endl;
    }
    tcolorer_.mpi_compute_range(localbodies_,range_,smoothinglength);
    return range_;
  }

  void update_iteration() 
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Destroy the previous tree
    if(tree_ !=  nullptr){
      delete tree_;
    }
    
  
    // Choose the smoothing length to be the biggest from everyone 
    double smoothinglength = 0;
    for(auto bi: localbodies_){
      if(smoothinglength < bi.second.getSmoothinglength()){
        smoothinglength = bi.second.getSmoothinglength();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&smoothinglength,1,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);

    if(rank==0){
      std::cout<<"H="<<smoothinglength<<std::endl;
    }

    // Then compute the range of the system 
    tcolorer_.mpi_compute_range(localbodies_,range_,smoothinglength);

    // Compute the keys 
    for(auto& bi:  localbodies_){
      bi.first = entity_key_t(range_,bi.second.coordinates());
    }

    // Distributed qsort and bodies exchange 
    tcolorer_.mpi_qsort(localbodies_,totalnbodies_);
    // Check for duplicates 
    assert(localbodies_.end() == std::unique(localbodies_.begin(),
          localbodies_.end(),[](const auto& left, const auto& right){
            return left.second.coordinates()==right.second.coordinates() &&
              left.first == right.first;
          }));
    
    // Generate the tree 
    tree_ = new tree_topology_t(range_[0],range_[1]);

    // Add my local bodies in my tree 
    // Clear the bodies_ vector 
    bodies_.clear();
    for(auto& bi:  localbodies_){
      auto nbi = tree_->make_entity(bi.second.getPosition(),&(bi.second),rank,
          bi.second.getMass());
      tree_->insert(nbi);
      bodies_.push_back(nbi);
    }

    // Exchnage usefull body_holder from my tree to other processes
    tcolorer_.mpi_branches_exchange(*tree_,localbodies_,rangeposproc_,
        range_,smoothinglength);

    // Compute and refresh the ghosts 
    tcolorer_.mpi_compute_ghosts(*tree_,smoothinglength,range_);
    //std::cout<<tree_->entities().size()<<std::endl;
    tcolorer_.mpi_refresh_ghosts(*tree_,range_);
  }

  void update_neighbors()
  {
    tcolorer_.mpi_refresh_ghosts(*tree_,range_);
  }

  void gravitation_fmm(){
    // Compute the COM
    tcolorer_.tree_traversal_com(*tree_);
    // Exchange the cells up to a depth 
    tcolorer_.mpi_exchange_cells(*tree_,maxmasscell_);
    // Compute the fmm interaction to the gathered cells 
    tcolorer_.mpi_compute_fmm(*tree_,macangle_);
    // Gather all the contributions and compute 
    tcolorer_.mpi_gather_cells(*tree_);
  }

  template<
    typename EF,
    typename... ARGS
  >
  void apply_in_smoothinglength(
      EF&& ef,
      ARGS&&... args)
  {
    for(auto& bi : bodies_){
      auto ents = tree_->find_in_radius(bi->getBody()->coordinates(),
          2*bi->getBody()->getSmoothinglength()/*+epsilon_*/);
      auto vecents = ents.to_vec();
      // Remove outside the h 
      //for(auto it = vecents.begin(); it != vecents.end() ; ){
        //if(flecsi::distance(bi->coordinates(),(*it)->coordinates())
        //      > 2.*bi->getBody()->getSmoothinglength()){
        //  it = vecents.erase(it);
        //}else{
        //  ++it;
        //}
        
      //}
      ef(bi,vecents,std::forward<ARGS>(args)...);
    } 
  }

  template<
    typename EF,
    typename... ARGS
  >
  void apply_all(
      EF&& ef,
      ARGS&&... args)
  {
    for(auto& bi: bodies_){
      ef(bi,std::forward<ARGS>(args)...);
    }
  }

private:
  int64_t totalnbodies_; 
  int64_t localnbodies_;
  double macangle_;
  double maxmasscell_;
  std::vector<std::pair<entity_key_t,body>> localbodies_;
  std::array<point_t,2> range_;
  std::vector<std::pair<point_t,point_t>> rangeposproc_;
  tree_colorer<T,D> tcolorer_;
  tree_topology_t* tree_;
  std::vector<body_holder*> bodies_;
  //double epsilon_ = 1.0;
};

#endif
