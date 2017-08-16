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
  }

  void read_bodies_txt(
      const char* filename)
  {
    io::inputDataTxtRange(localbodies_,localnbodies_,totalnbodies_,filename);
    totalmass_ = 0.;
    // Also compute the total mass 
    for(auto bi: localbodies_){
      totalmass_ += bi.second.getMass(); 
    }
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
    for(auto bi: localbodies_){
      if(smoothinglength_ < bi.second.getSmoothinglength()){
        smoothinglength_ = bi.second.getSmoothinglength();
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

    // Choose the smoothing length to be the biggest from everyone 
    smoothinglength_ = 0;
    for(auto bi: localbodies_){
      if(smoothinglength_ < bi.second.getSmoothinglength()){
        smoothinglength_ = bi.second.getSmoothinglength();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,&smoothinglength_,1,MPI_DOUBLE,MPI_MAX,
        MPI_COMM_WORLD);

    if(rank==0){
      std::cout<<"H="<<smoothinglength_<<std::endl;
    }
    tcolorer_.mpi_compute_range(localbodies_,range_,smoothinglength_);
    return range_;
  }

  void 
  update_iteration() 
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef DEBUG
    // Check for duplicate in the local bodies 
    auto tmp1 = localbodies_;
    int64_t colision = 0L;
    auto tmpit = std::unique(tmp1.begin(),tmp1.end(),
      [&rank](const auto& left, const auto& right){
        if(left.first == right.first){
          std::cout<<rank<<": "<<left.second<<" && "<<right.second<<std::endl;
          return true; 
        }
        return false; 
      });
    if(tmpit != tmp1.end()){
      std::cout<<rank<<": #colisions="<<std::distance(tmpit,tmp1.end())
      <<std::endl;
    }
#endif 

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

    // Compute the keys 
    for(auto& bi:  localbodies_){
      bi.first = entity_key_t(range_,bi.second.coordinates());
    }

#ifdef DEBUG
    double checkmassnt = 0.;
    for(auto bi: localbodies_){
      checkmassnt += bi.second.getMass(); 
    }
    MPI_Allreduce(MPI_IN_PLACE,&checkmassnt,1,
    MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    std::cout<<rank<<": "<< std::setprecision(10) <<
    checkmassnt<<" == "<<totalmass_<<" diff:"<<totalmass_-checkmassnt
    <<" min:"<<1.0e-2*minmass_<<std::endl<<std::flush;
    assert(fabs(checkmassnt-totalmass_) < 1.0e-2*minmass_); 
#endif   
 
    // Distributed qsort and bodies exchange 
    tcolorer_.mpi_qsort(localbodies_,totalnbodies_);
 
#ifdef DEBUG
    checkmassnt = 0.;
    for(auto bi: localbodies_){
      checkmassnt += bi.second.getMass(); 
    }
    MPI_Allreduce(MPI_IN_PLACE,&checkmassnt,1,
    MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    std::cout<<rank<<": "<< std::setprecision(10) <<
    checkmassnt<<" == "<<totalmass_<<" diff:"<<totalmass_-checkmassnt
    <<" min:"<<1.0e-2*minmass_<<std::endl<<std::flush;
    assert(fabs(checkmassnt-totalmass_) < 1.0e-2*minmass_); 
#endif   
 
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

#ifdef DEBUG
    double checkmass = 0;
    auto vect= tree_->entities().to_vec();
    for(auto v: vect){
      checkmass += v->getMass();
    } 
    std::cout<<rank<<": local="<<checkmass<<std::endl;
    MPI_Allreduce(MPI_IN_PLACE,&checkmass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    std::cout<<rank<<": av update "<< std::setprecision(10) <<
    checkmass<<" == "<<totalmass_<<" diff:"<<totalmass_-checkmass<<
    std::endl<<std::flush;
    assert(fabs(checkmass-totalmass_) < 1.0e-2*minmass_); 
#endif



    // Check the total number of bodies 
    int64_t checknparticles = bodies_.size();
    MPI_Allreduce(MPI_IN_PLACE,&checknparticles,1,MPI_INT64_T,
    MPI_SUM,MPI_COMM_WORLD); 
    assert(checknparticles==totalnbodies_);

    tree_->update_branches(2*smoothinglength_); 
    // Check the total mass of system 

#ifdef DEBUG
    checkmass = tree_->root()->getMass(); 
    std::cout<<rank<<": local="<<checkmass<<std::endl;
    MPI_Allreduce(MPI_IN_PLACE,&checkmass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    std::cout<<rank<<": ap update "<< std::setprecision(10) <<
    checkmass<<" == "<<totalmass_<<" diff:"<<totalmass_-checkmass<<
    std::endl<<std::flush;
    assert(fabs(checkmass-totalmass_) < 1.0e-2*minmass_); 
#endif

    // Exchnage usefull body_holder from my tree to other processes
    tcolorer_.mpi_branches_exchange(*tree_,localbodies_,rangeposproc_,
        range_,smoothinglength_);

    // Update the tree 
    tree_->update_branches(2*smoothinglength_); 
    //std::cout<<"TWO=="<<rank<<": "<<tree_->root()->getMass()<<std::endl;

    // Compute and refresh the ghosts 
    tcolorer_.mpi_compute_ghosts(*tree_,smoothinglength_,range_);
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
      for(size_t d=0; d<gdimension; ++d){
        assert(!std::isnan(bi->getBody()->coordinates()[d])); 
      }
      assert(bi->getBody()->getSmoothinglength() > 0.); 
      auto ents = tree_->find_in_radius_b(
        bi->getBody()->coordinates(),
        2*bi->getBody()->getSmoothinglength());
      auto vecents = ents.to_vec();
      if(vecents.size() == 0){
        std::cout<< "Particle:" << *(bi->getBody()) << std::endl 
        << "Holder:"<< *bi <<std::endl;
      }   
      assert(vecents.size()>0);

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
 //     if(bi->getBody()->getType() == 0){
        ef(bi,std::forward<ARGS>(args)...);
 //     }
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
