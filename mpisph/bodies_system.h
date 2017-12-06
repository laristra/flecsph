/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 
#ifndef _mpisph_body_system_h_
#define _mpisph_body_system_h_

#include <omp.h>

#include "tree_colorer.h"
#include "physics.h"
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
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    // Display the number of threads in DEBUG mode
    if(rank==0)
    {
      #pragma omp parallel 
      #pragma omp single 
      std::cout<<"OMP: "<<omp_get_num_threads()<<std::endl;
    }
  };

  /**
   * @brief      Destroys the object.
   */
  ~body_system(){
    if(tree_ != nullptr){
      delete tree_;
    }
  };

  // 
  
  /**
   * @brief      Max mass to stop the tree search during 
   * the gravitation computation with FMM
   *
   * @param[in]  maxmasscell  The maximum mass for the cells
   */
  void setMaxmasscell(
    double maxmasscell)
  {
    maxmasscell_ = maxmasscell;
  };
  
  /**
   * @brief      Sets the Multipole Acceptance Criterion for FMM
   *
   * @param[in]  macangle  Multipole Acceptance Criterion
   */
  void 
  setMacangle(
    double macangle)
  {
    macangle_ = macangle;
  };

  
  template<
    typename TL
  >
  TL
  get_attribute(
      const char * filename,
      const char * attributeName,
      TL default_value = TL(0))
  {
    TL value = TL{};
    if(typeid(TL)==typeid(double)){
      value = io::input_parameter_double(filename,attributeName);
    }else if(typeid(TL) == typeid(int)){
      value = io::input_parameter_int(filename,attributeName);
    }
    if(value == TL{}){
      value = default_value;
    }
    return value;
  }

  /**
   * @brief      Read the bodies from H5part file
   * Compute also the total to check for mass lost 
   *
   * @param[in]  filename        The filename
   * @param[in]  startiteration  The iteration from which load the data
   */
  void 
  read_bodies(
      const char * filename,
      int startiteration)
  {

    //io::init_reading(localbodies_,filename,totalnbodies_,localnbodies_,
    //    startiteration);
    //physics::read_data(filename,localbodies_,startiteration);
    io::inputDataHDF5(localbodies_,filename,totalnbodies_,localnbodies_,
        startiteration);
    
    #ifdef DEBUG
    minmass_ = 1.0e50;
    totalmass_ = 0.;
    // Also compute the total mass 
    for(auto bi: localbodies_){
      totalmass_ += bi.second.getMass(); 
      if(bi.second.getMass() < minmass_){
        minmass_ = bi.second.getMass(); 
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&minmass_,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE,&totalmass_,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
    #endif
  }

  /**
   * @brief      Write bodies to file in parallel 
   * Caution provide the file name prefix, h5part will be added
   * This is useful in case of multiple files output
   *
   * @param[in]  filename       The outut file prefix 
   * @param[in]  iter           The iteration of output
   * @param[in]  do_diff_files  Generate a file for each steps 
   */
  void 
  write_bodies(
      const char * filename, 
      int iter,
      bool do_diff_files = false)
  {
    io::outputDataHDF5(localbodies_,filename,iter,do_diff_files);
  }


  /**
   * @brief      Compute the largest smoothing length in the system 
   * This is really useful for particles with differents smoothing length
   *
   * @return     The largest smoothinglength of the system.
   */
  double 
  getSmoothinglength()
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Choose the smoothing length to be the biggest from everyone 
    smoothinglength_ = 0;
    for(auto bi: localbodies_){
      if(smoothinglength_ < bi.second.getSmoothinglength()){
        smoothinglength_ = bi.second.getSmoothinglength();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&smoothinglength_,1,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);

    #ifdef DEBUG
    #ifdef OUTPUT
    if(rank==0){
      std::cout<<"H="<<smoothinglength_<<std::endl;
    }
    #endif
    #endif 

    return smoothinglength_;

  }

  /**
   * @brief      Compute the range of the total particles 
   *
   * @return     The range.
   */
  std::array<point_t,2>& 
  getRange()
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Choose the smoothing length to be the biggest from everyone 
    smoothinglength_ = 0;
    for(auto bi: localbodies_){
      if(smoothinglength_ < bi.second.getSmoothinglength()){
        smoothinglength_ = bi.second.getSmoothinglength();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&smoothinglength_,1,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);

    #ifdef DEBUG
    #ifdef OUTPUT
    if(rank==0){
      std::cout<<"H="<<smoothinglength_<<std::endl;
    }
    #endif
    #endif
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
    for(auto& bi:  localbodies_){
      bi.first = entity_key_t(/*range_,*/bi.second.coordinates());
    }

 
    // Distributed qsort and bodies exchange 
    tcolorer_.mpi_qsort(localbodies_,totalnbodies_,neighbors_count_);
 
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
    std::cout<<rank<<" sub_entities="<<tree_->root()->sub_entities()<<std::endl<<std::flush; 
    //std::cout<<"TWO=="<<rank<<": "<<tree_->root()->getMass()<<std::endl;

#if COMPUTE_NEIGHBORS == 1 
    tcolorer_.mpi_compute_neighbors(
        *tree_,
        bodies_,
        neighbors_,
        neighbors_count_,
        smoothinglength_);
    // Compute and refresh the ghosts 
    tcolorer_.mpi_compute_ghosts(
        *tree_,
        bodies_,
        neighbors_,
        neighbors_count_,
        smoothinglength_/*,range_*/);
    tcolorer_.mpi_refresh_ghosts(*tree_/*,range_*/); 
#else 
// Compute and refresh the ghosts 
    tcolorer_.mpi_compute_ghosts(*tree_,bodies_,smoothinglength_/*,range_*/);
    tcolorer_.mpi_refresh_ghosts(*tree_/*,range_*/); 
#endif 
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

  #if COMPUTE_NEIGHBORS == 1 
  int64_t nelem = bodies_.size(); 
  #pragma omp parallel for 
    for(int64_t i=0; i<nelem; ++i){
      for(size_t d=0; d<gdimension; ++d){
        assert(!std::isnan(bodies_[i]->getBody()->coordinates()[d])); 
      }
      assert(bodies_[i]->getBody()->getSmoothinglength() > 0.); 
      assert(neighbors_count_[i+1]-neighbors_count_[i]>0);
      std::vector<body_holder*> nbs(neighbors_count_[i+1]-neighbors_count_[i]);
      int local = 0;
      for(int j=neighbors_count_[i]; j<neighbors_count_[i+1]; ++j)
      {
        nbs[local++]=neighbors_[j];
      }
      ef(bodies_[i],nbs,std::forward<ARGS>(args)...);
    } 
#else 
  int64_t nelem = bodies_.size(); 
  #pragma omp parallel for 
    for(int64_t i=0; i<nelem; ++i){
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
#endif
  }

  template<
    typename EF,
    typename... ARGS
  >
  void apply_all(
      EF&& ef,
      ARGS&&... args)
  {
    int64_t nelem = bodies_.size(); 
    #pragma omp parallel for 
    for(int64_t i=0; i<nelem; ++i){
        ef(bodies_[i],std::forward<ARGS>(args)...);
    }
  }

  template<
    typename EF,
    typename... ARGS
  >
  void get_all(
    EF&& ef, 
    ARGS&&... args)
  {
    ef(bodies_,std::forward<ARGS>(args)...);
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
  std::vector<std::array<point_t,2>> rangeposproc_;
  tree_colorer<T,D> tcolorer_;
  tree_topology_t* tree_;
  std::vector<body_holder*> bodies_;
  double smoothinglength_;
  double totalmass_;
  double minmass_;
  
  std::vector<int64_t> neighbors_count_;
  std::vector<body_holder*> neighbors_; 
  //double epsilon_ = 1.0;
};

#endif
