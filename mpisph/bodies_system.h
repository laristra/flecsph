/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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

#include "fmm.h"
#include "io.h"
#include "params.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <omp.h>
#include <typeinfo>

#include "psort.h"

#define DEBUG_TREE

using namespace mpi_utils;

/**
 * @brief      The bodies/particles system.
 * This is a wrapper for a simpler use from users.
 *
 * @tparam     T     The type of data, usualy double
 * @tparam     D     The dimension of the current simulation
 */
template<typename T, size_t D>
class body_system
{

  using point_t = flecsi::space_vector_u<T, D>;

public:
  /**
   * @brief      Constructs the object.
   */
  body_system()
    : totalnbodies_(0L), localnbodies_(0L), macangle_(0.0),
      maxmasscell_(1.0e-40), tree_{} {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Display the number of threads in DEBUG mode

    if(param::sph_variable_h) {
      log_one(warn) << "Variable smoothing length ENABLE" << std::endl;
    }
  };

  /**
   * @brief      Destroys the object.
   */
  ~body_system(){};

  /**
   * @brief      Sets the Multipole Acceptance Criterion for FMM
   *
   * @param[in]  macangle  Multipole Acceptance Criterion
   */
  void setMacangle(double macangle) {
    macangle_ = macangle;
  };

  /**
   * @brief      Read the bodies from H5part file Compute also the total to
   *             check for mass lost
   *
   * @param[in]  input_prefix    Input filename without format extension or
   *                             step number (i.e. sim_00000.h5part -> "sim")
   * @param[in]  output_prefix   Output filename prefix
   * @param[in]  startiteration  The iteration from which load the data
   */
  void read_bodies(const char * input_prefix,
    const char * output_prefix,
    const int startiteration) {

    io::inputDataHDF5(tree_.entities(), input_prefix, output_prefix,
      totalnbodies_, localnbodies_, startiteration);
  }

  /**
   * @brief      Write bodies to file in parallel Caution provide the
   * file name
   *             prefix, h5part will be added This is useful in case
   *             of multiple
   *             files output
   *
   * @param[in]  output_prefix  The output file prefix
   * @param[in]  iter           The iteration of output
   * @param[in]  do_diff_files  Generate a file for each steps
   */
  void write_bodies(const char * output_prefix, int iter, double totaltime) {
    io::outputDataHDF5(tree_.entities(), output_prefix, iter, totaltime);
  }

  /**
   * @brief      Compute the largest smoothing length in the system This is
   *             really useful for particles with differents smoothing length
   *
   * @return     The largest smoothinglength of the system.
   */
  double getSmoothinglength() {
    // Choose the smoothing length to be the biggest from everyone
    double smoothinglength = 0;
    for(size_t i = 0; i < tree_.entities().size(); ++i) {
      if(smoothinglength < tree_.entity(i).radius()) {
        smoothinglength = tree_.entity(i).radius();
      }
    }

    MPI_Allreduce(
      MPI_IN_PLACE, &smoothinglength, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return smoothinglength;
  }

  /**
   * @brief      Compute the range of thw whole particle system
   *
   * @return     The range.
   */
  std::array<point_t, 2> & getRange() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    mpi_compute_range(tree_.entities(), range_);
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
  void update_iteration() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Clean the whole tree structure
    tree_.clean();

    if(param::periodic_boundary_x || param::periodic_boundary_y ||
       param::periodic_boundary_z) {
      boundary::pboundary_clean(tree_.entities());
      // Choose the smoothing length to be the biggest from everyone
      double smoothinglength = getSmoothinglength();
      boundary::pboundary_generate(tree_.entities(), 2.5 * smoothinglength);
      localnbodies_ = tree_.entities().size();
      MPI_Allreduce(&localnbodies_, &totalnbodies_, 1, MPI_INT64_T, MPI_SUM,
        MPI_COMM_WORLD);
    }

    log_one(trace) << "#particles: " << totalnbodies_ << std::endl;
    // Then compute the range of the system
    mpi_compute_range(tree_.entities(), range_);
    if(range_[0] == range_[1]) {
      std::cerr << "Range are equals: " << range_[0] << " == " << range_[1]
                << std::endl;
      assert(range_[0] != range_[1]);
    }
    log_one(trace) << "Range=" << range_[0] << std::endl;
    log_one(trace) << "      " << range_[1] << std::endl;
    // Generate the tree based on the range
    tree_.set_range(range_);
    // Compute the keys
    tree_.compute_keys();

    // Distributed sort
    log_one(trace) << "QSort (" << size << ")" << std::endl;
    double timer = omp_get_wtime();

    int dist[size];
    dist[rank] = tree_.entities().size();

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, dist, 1, MPI_INT, MPI_COMM_WORLD);

    psort::psort(
      tree_.entities(),
      [](auto & left, auto & right) {
        if(left.key() < right.key()) {
          return true;
        }
        if(left.key() == right.key()) {
          return left.id() < right.id();
        }
        return false;
      },
      dist);
    log_one(trace) << "QSort.done: ppp=" << tree_.entities().size() << "+-1 "
                   << omp_get_wtime() - timer << "s" << std::endl;

#ifdef DEBUG_TREE
    std::vector<int> totalprocbodies;
    totalprocbodies.resize(size);
    int mybodies = tree_.entities().size();
    // Share the final array size of everybody
    MPI_Allgather(
      &mybodies, 1, MPI_INT, &totalprocbodies[0], 1, MPI_INT, MPI_COMM_WORLD);
    int min = *std::min_element(totalprocbodies.begin(), totalprocbodies.end());
    int max = *std::max_element(totalprocbodies.begin(), totalprocbodies.end());
    int total = std::accumulate(totalprocbodies.begin(), totalprocbodies.end(), 0); 
    assert(total == totalnbodies_); 
    assert(max - min <= 1);
#endif // DEBUG_TREE

    tree_.build_tree(physics::compute_cofm);
    log_one(trace) << "#particles: " << totalnbodies_ << std::endl;

    localnbodies_ = tree_.entities().size();
    log_one(trace) << tree_ << std::endl;
  }

  void mpi_compute_range(const std::vector<body> & bodies,
    std::array<point_t, 2> & range) {
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Compute the local range
    range_t lrange;

    lrange[1] = bodies.back().coordinates();
    lrange[0] = bodies.back().coordinates();

    for(size_t i = 0; i < bodies.size(); ++i) {
      for(size_t d = 0; d < gdimension; ++d) {
        if(bodies[i].coordinates()[d] + bodies[i].radius() > lrange[1][d])
          lrange[1][d] = bodies[i].coordinates()[d] + bodies[i].radius();
        if(bodies[i].coordinates()[d] - bodies[i].radius() < lrange[0][d])
          lrange[0][d] = bodies[i].coordinates()[d] - bodies[i].radius();
      } // for
    } // for

    double max[gdimension];
    double min[gdimension];
    for(size_t i = 0; i < gdimension; ++i) {
      max[i] = lrange[1][i];
      min[i] = lrange[0][i];
    } // for

    MPI_Allreduce(
      MPI_IN_PLACE, max, gdimension, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(
      MPI_IN_PLACE, min, gdimension, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    for(size_t d = 0; d < gdimension; ++d) {
      range[0][d] = min[d];
      range[1][d] = max[d];
    } // for
  }

  /**
   * Reset the ghosts of the tree to start in the next tree traversal
   */
  void reset_ghosts() {
    tree_.reset_ghosts(physics::compute_cofm);
  }

  /**
   * @brief      Compute the gravition interction between all the particles
   * @details    The function is based on Fast Multipole Method. The functions
   *             are defined in the file tree_fmm.h
   */
  void gravitation_fmm() {
    assert (gdimension == 3);
    if constexpr (gdimension == 3) { 
      using namespace fmm;
      tree_.traversal_fmm(macangle_, taylor_c2c, taylor_p2c, fmm_p2p, fmm_c2p);
    }
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
  template<typename EF, typename... ARGS>
  void apply_in_smoothinglength(EF && ef, ARGS &&... args) {
    tree_.traversal_sph(ef, std::forward<ARGS>(args)...);
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
  template<typename EF, typename... ARGS>
  void apply_all(EF && ef, ARGS &&... args) {
    int64_t nelem = tree_.entities().size();
    for(int64_t i = 0; i < nelem; ++i) {
      ef(tree_.entities()[i], std::forward<ARGS>(args)...);
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
  template<typename EF, typename... ARGS>
  void get_all(EF && ef, ARGS &&... args) {
    ef(tree_.entities(), std::forward<ARGS>(args)...);
  }

  /**
   * @brief      Gets a vector of the local bodies of this process.
   *
   * @return     The localbodies.
   */
  std::vector<body> & getLocalbodies() {
    return tree_.entities();
  };

  size_t nbodies() {
    return tree_.entities().size();
  }

  /**
   * @ brief return the number of local bodies
   */
  int64_t getNLocalBodies() {
    return localnbodies_;
  }

  int64_t getNBodies() {
    return totalnbodies_;
  }

  tree_topology_t * tree() {
    return &tree_;
  }

private:
  int64_t totalnbodies_; // Total number of local particles
  int64_t localnbodies_; // Local number of particles
  double macangle_; // Macangle for FMM
  double maxmasscell_; // Mass criterion for FMM
  range_t range_;
  tree_topology_t tree_; // The particle tree data structure
  double epsilon_ = 0.;

  const int refresh_tree = 0;
  int current_refresh = refresh_tree;
};

#endif
