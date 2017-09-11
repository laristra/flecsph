/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 
#ifndef _mpisph_body_system_h_
#define _mpisph_body_system_h_

#include "tree_colorer.h"
#include "io.h"
#include "omp.h"

template<
  typename T,
  size_t D
  >
class body_system{

using point_t = flecsi::point<T,D>;

public:
  body_system():totalnbodies_(0L),localnbodies_(0L),macangle_(0.0),
  maxmasscell_(1.0e-40),tree_(nullptr)
  {
#ifdef DEBUG
#pragma omp parallel  
#pragma omp single 
    std::cout<<"OMP: "<<omp_get_num_threads()<<std::endl;
#endif
  };
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
    minmass_ = 1.0e50;
    totalmass_ = 0.;

    // Also compute the total mass 
#pragma omp parallel for reduction(+:totalmass_) reduction(min: minmass_)
    for(int i = 0; i < localnbodies_; ++i){
      totalmass_ += localbodies_[i].second.getMass(); 
      if(localbodies_[i].second.getMass() < minmass_){
        minmass_ = localbodies_[i].second.getMass(); 
      }
    }

    MPI_Allreduce(MPI_IN_PLACE,&minmass_,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE,&totalmass_,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
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
    smoothinglength_ = 0;

#pragma omp parallel for reduction(max: smoothinglength_)
    for(int i=0; i < localnbodies_; ++i){
      if(smoothinglength_ < localbodies_[i].second.getSmoothinglength()){
        smoothinglength_ = localbodies_[i].second.getSmoothinglength();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&smoothinglength_,1,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);

    if(rank==0){
      std::cout<<"H="<<smoothinglength_<<std::endl;
    }
    return smoothinglength_;

  }

  std::array<point_t,2>& 
  getRange()
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    smoothinglength_ = getSmoothinglength(); 
    
    tcolorer_.mpi_compute_range(localbodies_,range_,smoothinglength_);
    return range_;
  }

  void 
  update_iteration() 
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Destroy the previous tree
    if(tree_ !=  nullptr){
      delete tree_;
    }
     
    // Choose the smoothing length to be the biggest from everyone 
    smoothinglength_ = getSmoothinglength();

    if(rank==0){
      std::cout<<"H="<<smoothinglength_<<std::endl;
    }

    // Then compute the range of the system 
    tcolorer_.mpi_compute_range(localbodies_,range_,smoothinglength_);


    // Setup the keys range 
    entity_key_t::set_range(range_); 
    
    // Compute the keys 
#pragma omp parallel for 
    for(int i=0; i<localnbodies_; ++i){
      localbodies_[i].first = entity_key_t(
          localbodies_[i].second.coordinates());
    }

    // Distributed qsort and bodies exchange 
    tcolorer_.mpi_qsort(localbodies_,totalnbodies_);
 
    // Generate the tree 
    tree_ = new tree_topology_t(range_[0],range_[1]);

    // Add my local bodies in my tree 
    // Clear the bodies_ vector 
    bodies_.clear();
    for(auto& bi:  localbodies_){
      auto nbi = tree_->make_entity(bi.second.getPosition(),&(bi.second),rank,
          bi.second.getMass(),bi.second.getId());
      tree_->insert(nbi);
      bodies_.push_back(nbi);
    }

    // Check the total number of bodies 
    int64_t checknparticles = bodies_.size();
    MPI_Allreduce(MPI_IN_PLACE,&checknparticles,1,MPI_INT64_T,
    MPI_SUM,MPI_COMM_WORLD); 
    assert(checknparticles==totalnbodies_);

    tree_->update_branches(2*smoothinglength_); 
    // Check the total mass of system 

    // Exchnage usefull body_holder from my tree to other processes
    tcolorer_.mpi_branches_exchange(*tree_,localbodies_,rangeposproc_,
        range_,smoothinglength_);

    // Update the tree 
    tree_->update_branches(2*smoothinglength_); 

    // Compute and refresh the ghosts 
    tcolorer_.mpi_compute_ghosts(*tree_,smoothinglength_/*,range_*/);
    tcolorer_.mpi_refresh_ghosts(*tree_/*,range_*/); 
  }

  void update_neighbors()
  {
    tcolorer_.mpi_refresh_ghosts(*tree_/*,range_*/);
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

#pragma omp parallel for 
    for(int i=0; i<localnbodies_; ++i){

      for(size_t d=0; d<gdimension; ++d){
        assert(!std::isnan(bodies_[i]->getBody()->coordinates()[d])); 
      }
      assert(bodies_[i]->getBody()->getSmoothinglength() > 0.); 
      auto ents = tree_->find_in_radius_b(
        bodies_[i]->getBody()->coordinates(),
        2*bodies_[i]->getBody()->getSmoothinglength());
      auto vecents = ents.to_vec();
      if(vecents.size() == 0){
        std::cout<< "Particle:" << *(bodies_[i]->getBody()) << std::endl 
        << "Holder:"<< *bodies_[i] <<std::endl;
      }   
      assert(vecents.size()>0);

      ef(bodies_[i],vecents,std::forward<ARGS>(args)...);
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
#pragma omp parallel for 
    for(int i=0; i<localnbodies_; ++i){
      ef(bodies_[i],std::forward<ARGS>(args)...);
    }
  }


  std::vector<std::pair<entity_key_t,body>>& 
    getLocalbodies(){
    return localbodies_;
  };

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
  double smoothinglength_;
  double totalmass_;
  double minmass_;
  //double epsilon_ = 1.0;
};

#endif
