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
#include "tree_colorer.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <omp.h>
#include <typeinfo>

using namespace mpi_utils;

//#define OUTPUT_TREE_GRAPH 1

/**
 * @brief      The bodies/particles system.
 * This is a wrapper for a simpler use from users.
 *
 * @tparam     T     The type of data, usualy double
 * @tparam     D     The dimension of the current simulation
 */
template <typename T, size_t D> class body_system {

  using point_t = flecsi::point_u<T, D>;

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

#pragma omp parallel
#pragma omp master
    clog_one(warn) << "USING OMP THREADS: " << omp_get_num_threads()
                   << std::endl;

    if (param::sph_variable_h) {
      clog_one(warn) << "Variable smoothing length ENABLE" << std::endl;
    }
  };

  /**
   * @brief      Destroys the object.
   */
  ~body_system(){};

  /**
   * @brief      Max mass to stop the tree search during
   * the gravitation computation with FMM
   *
   * @param[in]  maxmasscell  The maximum mass for the cells
   */
  void setMaxmasscell(double maxmasscell) { maxmasscell_ = maxmasscell; };

  /**
   * @brief      Sets the Multipole Acceptance Criterion for FMM
   *
   * @param[in]  macangle  Multipole Acceptance Criterion
   */
  void setMacangle(double macangle) { macangle_ = macangle; };

  /**
   * @brief      Read the bodies from H5part file Compute also the total to
   *             check for mass lost
   *
   * @param[in]  input_prefix    Input filename without format extension or
   *                             step number (i.e. sim_00000.h5part -> "sim")
   * @param[in]  output_prefix   Output filename prefix
   * @param[in]  startiteration  The iteration from which load the data
   */
  void read_bodies(const char *input_prefix, const char *output_prefix,
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
  void write_bodies(const char *output_prefix, int iter, double totaltime) {
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
#pragma omp parallel for reduction(max : smoothinglength)
    for (size_t i = 0; i < tree_.entities().size(); ++i) {
      if (smoothinglength < tree_.entity(i).radius()) {
        smoothinglength = tree_.entity(i).radius();
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &smoothinglength, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    return smoothinglength;
  }

  /**
   * @brief      Compute the range of thw whole particle system
   *
   * @return     The range.
   */
  std::array<point_t, 2> &getRange() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    tcolorer_.mpi_compute_range(tree_.entities(), range_);
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
    if (current_refresh == refresh_tree) {
      clog_one(trace) << "Reset tree" << std::endl;
      tree_.clean();
    } else {
      clog_one(trace) << "Reset Ghosts" << std::endl;
      tree_.reset_ghosts(false);
    }

    if (param::periodic_boundary_x || param::periodic_boundary_y ||
        param::periodic_boundary_z) {
      boundary::pboundary_clean(tree_.entities());
      // Choose the smoothing length to be the biggest from everyone
      double smoothinglength = getSmoothinglength();
      boundary::pboundary_generate(tree_.entities(), 2.5 * smoothinglength);
      localnbodies_ = tree_.entities().size();
      MPI_Allreduce(&localnbodies_, &totalnbodies_, 1, MPI_INT64_T, MPI_SUM,
                    MPI_COMM_WORLD);
    }

    clog_one(trace) << "#particles: " << totalnbodies_ << std::endl;

    if (current_refresh == refresh_tree) {
      clog_one(trace) << "Exchange everythings" << std::endl;
      // Then compute the range of the system
      tcolorer_.mpi_compute_range(tree_.entities(), range_);
      assert(range_[0] != range_[1]);
      clog_one(trace) << "Range=" << range_[0] << std::endl;
      clog_one(trace) << "      " << range_[1] << std::endl;
      // Generate the tree based on the range
      tree_.set_range(range_);
      // Compute the keys
      tree_.compute_keys();
      // Distributed sample sort
      tcolorer_.mpi_qsort(tree_.entities(), totalnbodies_);
    }

#ifdef OUTPUT_TREE_INFO
    clog_one(trace) << "Construction of the tree";
#endif

    if (current_refresh == refresh_tree) {
// Sort the bodies
#ifdef BOOST_PARALLEL
      boost::sort::block_indirect_sort(
#else
      std::sort(
#endif
          tree_.entities().begin(), tree_.entities().end(),
          [](auto &left, auto &right) {
            if (left.key() < right.key()) {
              return true;
            }
            if (left.key() == right.key()) {
              return left.id() < right.id();
            }
            return false;
          }); // sort
    }

    if (current_refresh == refresh_tree) {
      // Add my local bodies in my tree
      // Clear the bodies_ vector
      for (auto &bi : tree_.entities()) {
        bi.set_owner(rank);
        auto id = tree_.make_entity(bi.key(), bi.coordinates(), &(bi), rank,
                                    bi.mass(), bi.id(), bi.radius());
        tree_.insert(id);
        auto nbi = tree_.get(id);
        assert(nbi->global_id() == bi.id());
        assert(nbi->entity_ptr() != nullptr);
        assert(nbi->is_local());
      }
      localnbodies_ = tree_.entities().size();
    }

#ifdef OUTPUT_TREE_INFO
    clog_one(trace) << ".done" << std::endl;
#endif
#ifdef OUTPUT_TREE_GRAPH
    tree_.mpi_tree_traversal_graphviz(0);
#endif

    if (!(param::periodic_boundary_x || param::periodic_boundary_y ||
          param::periodic_boundary_z)) {
#ifdef DEBUG
      // Check the total number of bodies
      int64_t checknparticles = tree_.tree_entities().size();
      MPI_Allreduce(MPI_IN_PLACE, &checknparticles, 1, MPI_INT64_T, MPI_SUM,
                    MPI_COMM_WORLD);
      assert(checknparticles == totalnbodies_);
#endif
    }
    // Add edge bodies from my direct neighbor
    tree_.share_edge();
    tree_.cofm(tree_.root(), epsilon_, false);
#ifdef OUTPUT_TREE_GRAPH
    tree_.mpi_tree_traversal_graphviz(1);
#endif

#ifdef OUTPUT_TREE_INFO
    {
      clog_one(trace) << "Computing branches" << std::endl;
      std::ostringstream oss;
      std::vector<int> nentities(size);
      int lentities = tree_.root()->sub_entities();
      // Get on 0
      MPI_Gather(&lentities, 1, MPI_INT, &nentities[0], 1, MPI_INT, 0,
                 MPI_COMM_WORLD);

      oss << rank << " sub_entities before=";
      for (auto v : nentities) {
        oss << v << ";";
      }
      oss << std::endl;
      clog_one(trace) << oss.str() << std::flush;

      oss.str("");
      oss.clear();
      clog_one(trace) << tree_ << std::endl;
    }
#endif

    // Exchnage usefull body_holder from my tree to other processes
    if (param::enable_fmm) {
      tcolorer_.mpi_branches_exchange_all_leaves(tree_, tree_.entities());
    } else {
      tcolorer_.mpi_branches_exchange(tree_, tree_.entities());
    }

    tree_.cofm(tree_.root(), epsilon_, false);
#ifdef OUTPUT_TREE_GRAPH
    tree_.mpi_tree_traversal_graphviz(2);
#endif

#ifdef OUTPUT_TREE_INFO
    {
      std::ostringstream oss;
      std::vector<int> nentities(size);
      int lentities = tree_.root()->sub_entities();
      // Get on 0
      MPI_Gather(&lentities, 1, MPI_INT, &nentities[0], 1, MPI_INT, 0,
                 MPI_COMM_WORLD);
      if (rank == 0) {
        oss << rank << " sub_entities after=";
        for (auto v : nentities) {
          oss << v << ";";
        }
        oss << std::endl;
        clog_one(trace) << oss.str() << std::flush;
      }
    }
#endif

#ifdef OUTPUT_TREE_INFO
    // Tree informations
    clog_one(trace) << tree_ << std::endl;
#endif
    if (current_refresh == 0) {
      current_refresh = refresh_tree;
    } else {
      --current_refresh;
    }
  }

  /**
   * Reset the ghosts of the tree to start in the next tree traversal
   */
  void reset_ghosts() { tree_.reset_ghosts(); }

  /**
   * @brief      Compute the gravition interction between all the particles
   * @details    The function is based on Fast Multipole Method. The functions
   *             are defined in the file tree_fmm.h
   */
  void gravitation_fmm() {
    tree_.traversal_fmm(tree_.root(), maxmasscell_, macangle_,
                        fmm::gravitation_fc, fmm::gravitation_dfcdr,
                        fmm::gravitation_dfcdrdr, fmm::interation_c2p);
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
  template <typename EF, typename... ARGS>
  void apply_in_smoothinglength(EF &&ef, ARGS &&... args) {
    tree_.traversal_sph(tree_.root(), ef, std::forward<ARGS>(args)...);
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
  template <typename EF, typename... ARGS>
  void apply_all(EF &&ef, ARGS &&... args) {
    int64_t nelem = tree_.entities().size();
#pragma omp parallel for
    for (int64_t i = 0; i < nelem; ++i) {
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
  template <typename EF, typename... ARGS>
  void get_all(EF &&ef, ARGS &&... args) {
    ef(tree_.entities(), std::forward<ARGS>(args)...);
  }

  /**
   * @brief      Test function using the n^2 algorithm testing
   *
   * @param[in]  <unnamed>  The function to apply
   * @param[in]  <unnamed>  The arguments of the function
   *
   * @tparam     EF         The function to apply
   * @tparam     ARGS       The arguments of the function
   */
  template <typename EF, typename... ARGS>
  void apply_square(EF &&ef, ARGS &&... args) {
    int64_t nelem = tree_.tree_entities().size();
#pragma omp parallel for
    for (int64_t i = 0; i < nelem; ++i) {
      ef(tree_.get(i), tree_.tree_entities(), std::forward<ARGS>(args)...);
    }
  }

  /**
   * @brief      Gets a vector of the local bodies of this process.
   *
   * @return     The localbodies.
   */
  std::vector<body> &getLocalbodies() { return tree_.entities(); };

  /**
   * @ brief return the number of local bodies
   */
  int64_t getNLocalBodies() { return localnbodies_; }

  int64_t getNBodies() { return totalnbodies_; }

  tree_topology_t *tree() { return &tree_; }

private:
  int64_t totalnbodies_; // Total number of local particles
  int64_t localnbodies_; // Local number of particles
  double macangle_;      // Macangle for FMM
  double maxmasscell_;   // Mass criterion for FMM
  range_t range_;
  std::vector<range_t> rangeposproc_;
  tree_colorer<T, D> tcolorer_;
  tree_topology_t tree_; // The particle tree data structure
  double epsilon_ = 0.;

  const int refresh_tree = 0;
  int current_refresh = refresh_tree;
};

#endif
