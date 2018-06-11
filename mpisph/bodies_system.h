/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

/**
 * @file bodies_system.h
 * @author Julien Loiseau 
 * @brief Class and function to handle the system of bodies/particles. 
 * Contain the function for user, hidding the IO/distribution and tree search.
 */
 
#ifndef _mpisph_body_system_h_
#define _mpisph_body_system_h_

#include "tree_colorer.h"
#include "tree_fmm.h"
#include "physics.h"
#include "io.h"
#include "utils.h"

#include <omp.h>
#include <iostream>
#include <fstream>

using namespace mpi_utils;

/**
 * @brief      The bodies/particles system. 
 * This is a wrapper for a simpler use from users. 
 *
 * @tparam     T     The type of data, usualy double
 * @tparam     D     The dimension of the current simulation
 */
template<
  typename T,
  size_t D
  >
class body_system{

using point_t = flecsi::point<T,D>;

public:

  /**
   * @brief      Constructs the object.
   */
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

  /**
   * @brief      Get the value of an attribute from an HDF5 file
   * @details    \TODO add more types of data 
   *
   * @param      filename       The file from which we get the attribute
   * @param      attributeName  The attribute name in char
   * @param      default_value  The default type 
   *
   * @tparam     TL             The attribute value 
   *
   * @return     The value of the attribute 
   */
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
   * @brief      Read the bodies from H5part file Compute also the total to
   *             check for mass lost
   *
   * @param[in]  filename        The filename
   * @param[in]  startiteration  The iteration from which load the data
   */
  void 
  read_bodies(
      const char * filename,
      int startiteration)
  {

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
   * @brief      Write bodies to file in parallel Caution provide the file name
   *             prefix, h5part will be added This is useful in case of multiple
   *             files output
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
   * @brief      Compute the largest smoothing length in the system This is
   *             really useful for particles with differents smoothing length
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
   * @brief      Compute the range of thw whole particle system 
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

  /**
   * @brief      Generate and share the particle for this iteration
   * @details    This part if decomposed with:
   *    - Compute and prepare the tree for this iteration
   *    - Compute the Max smoothing length
   *    - Compute the range of the system using the smoothinglength
   *    - Cmopute the keys
   *    - Distributed qsort and sharing
   *    - Generate and feed the tree
   *    - Exchange branches for smoothing length
   *    - Compute and exchange ghosts in real smoothing length
   */
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
    localnbodies_ = localbodies_.size();

    // Check the total number of bodies 
    int64_t checknparticles = bodies_.size();
    MPI_Allreduce(MPI_IN_PLACE,&checknparticles,1,MPI_INT64_T,
    MPI_SUM,MPI_COMM_WORLD); 
    assert(checknparticles==totalnbodies_);

    tree_->update_branches(2*smoothinglength_); 

#ifdef DEBUG
    std::vector<int> nentities(size);
    int lentities = tree_->root()->sub_entities();
    // Get on 0 
    MPI_Gather(
      &lentities,
      1,
      MPI_INT,
      &nentities[0],
      1,
      MPI_INT,
      0,
      MPI_COMM_WORLD
      );

    if(rank == 0){
      std::cout<<rank<<" sub_entities before="; 
      for(auto v: nentities){
        std::cout<<v<<";";
      }
      std::cout<<std::endl;
    }
#endif

    // Exchnage usefull body_holder from my tree to other processes
    tcolorer_.mpi_branches_exchange(*tree_,localbodies_,rangeposproc_,
        range_,smoothinglength_);

    // Update the tree 
    tree_->update_branches(2*smoothinglength_);

#ifdef DEBUG
    lentities = tree_->root()->sub_entities();
    // Get on 0 
    MPI_Gather(
      &lentities,
      1,
      MPI_INT,
      &nentities[0],
      1,
      MPI_INT,
      0,
      MPI_COMM_WORLD
      );

    if(rank == 0){
      std::cout<<rank<<" sub_entities after="; 
      for(auto v: nentities){
        std::cout<<v<<";";
      }
      std::cout<<std::endl;
    }
#endif
    
    tcolorer_.mpi_compute_ghosts(*tree_,bodies_,smoothinglength_/*,range_*/);
    tcolorer_.mpi_refresh_ghosts(*tree_/*,range_*/); 

  }

  /**
   * @brief      Update the neighbors that have beem compute in update_iteration
   * This function use buffer pre-computed to update the data faster. 
   */
  void update_neighbors()
  {
    tcolorer_.mpi_refresh_ghosts(*tree_/*,range_*/);
  }

  /**
   * @brief      Compute the gravition interction between all the particles
   * @details    The function is based on Fast Multipole Method. The functions
   *             are defined in the file tree_fmm.h
   */
  void 
  gravitation_fmm()
  {
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank == 0){
      std::cout<<"FMM: mmass="<<maxmasscell_<<" angle="<<macangle_<<std::endl;
    }

    // Just consider the local particles in the tree for FMM 
    tree_->update_branches_local(smoothinglength_);
    assert((int64_t)tree_->root()->sub_entities() == localnbodies_);

    tfmm_.mpi_exchange_cells(*tree_,maxmasscell_);
    tfmm_.mpi_compute_fmm(*tree_,macangle_,0);
    tfmm_.mpi_gather_cells(*tree_,macangle_,totalnbodies_);
    
    // Reset the tree to normal before leaving
    tree_->update_branches(2*smoothinglength_);
  }

  /**
   * @brief      Apply the function EF with ARGS in the smoothing length of all
   *             the lcoal particles. This function need a previous call to
   *             update_iteration and update_neighbors for the remote particles'
   *             data.
   *
   * @param[in]  ef    The function to apply in the smoothing length
   * @param[in]  args  Arguments of the physics function applied in the 
   *                   smoothing length
   *
   * @tparam     EF         The function to apply in the smoothing length
   * @tparam     ARGS       Arguments of the physics function applied in the
   *                        smoothing length
   */
  template<
    typename EF,
    typename... ARGS
  >
  void apply_in_smoothinglength(
      EF&& ef,
      ARGS&&... args)
  {
    int64_t ncritical = 32; 
    tree_->apply_sub_cells(
        tree_->root(),
        bodies_[0]->getBody()->getSmoothinglength()*2.,
        0.,
        ncritical,
        ef,
        std::forward<ARGS>(args)...); 
  }

  /**
   * @brief      Apply a function to all the particles. 
   *
   * @param[in]  <unnamed>  { parameter_description }
   * @param[in]  <unnamed>  { parameter_description }
   *
   * @tparam     EF         The function to apply to all particles
   * @tparam     ARGS       Arguments of the function for all particles
   */
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

  /**
   * @brief      Apply a function on the vector of local bodies 
   *
   * @param[in]  <unnamed>  { parameter_description }
   * @param[in]  <unnamed>  { parameter_description }
   *
   * @tparam     EF         The function to apply to the vector 
   * @tparam     ARGS       Arguments of the function to apply to the vector
   */
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

  /**
   * @brief      Gets a vector of the local bodies of this process.
   *
   * @return     The localbodies.
   */
  std::vector<std::pair<entity_key_t,body>>& 
    getLocalbodies(
      )
  {
    return localbodies_;
  };

private:
  int64_t totalnbodies_;        // Total number of local particles
  int64_t localnbodies_;        // Local number of particles
  double macangle_;             // Macangle for FMM
  double maxmasscell_;          // Mass criterion for FMM
  std::vector<std::pair<entity_key_t,body>> localbodies_;
  std::array<point_t,2> range_;
  std::vector<std::array<point_t,2>> rangeposproc_;
  tree_colorer<T,D> tcolorer_;
  tree_fmm<T,D> tfmm_;        // tree_fmm.h function for FMM 
  tree_topology_t* tree_;     // The particle tree data structure
  std::vector<body_holder*> bodies_;
  double smoothinglength_;    // Keep track of the biggest smoothing length 
  double totalmass_;          // Check the total mass of the system 
  double minmass_;            // Check the minimal mass of the system
  
  std::vector<int64_t> neighbors_count_;
  std::vector<body_holder*> neighbors_; 
  //double epsilon_ = 1.0;
};

#endif
