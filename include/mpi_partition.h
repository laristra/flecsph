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
 * @file mpi_partition.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Function needed for MPI distribution of the bodies 
 */

#ifndef mpi_partition_h
#define mpi_partition_h

#include <vector>
#include <numeric>
#include <iostream>

#include "tree.h"

std::ostream&
operator<<(
    std::ostream& ostr,
    const entity_key_t& id
);

struct mpi_cell{
  point_t position;
  point_t fc;
  double dfcdr[9];
  double dfcdrdr[27];
  point_t bmax;
  point_t bmin;
  //double radius;
  branch_id_t id;

  mpi_cell(point_t position_,point_t bmin_,point_t bmax_, branch_id_t id_){
    position = position_;  
    id = id_;
    bmax = bmax_;
    bmin = bmin_;
  }
  mpi_cell(){};
};

struct mpi_ghosts_t{
  std::vector<body> sendbodies;
  std::vector<int> nsendholders;
  std::vector<int> nsendoffsets;
  std::vector<body> recvbodies;
  std::vector<int> nrecvholders;
  std::vector<int> nrecvoffsets;
  std::vector<body_holder*> totalrecvholders;
  std::vector<std::set<body_holder*>> sendholders;
};

void
mpi_gather_ghosts_com(
    tree_topology_t&,
    std::vector<mpi_cell>&,
    std::vector<int>&,
    std::array<point_t,2>&);

void 
mpi_gather_cells(
    tree_topology_t&,
    std::vector<mpi_cell> &,
    std::vector<int>&);

void 
mpi_compute_fmm(
    tree_topology_t&,
    std::vector<mpi_cell>&,
    double);

void 
mpi_exchange_cells(
    tree_topology_t&, 
    std::vector<mpi_cell>&,
    std::vector<int>&,
    double);

void 
mpi_compute_range(
  std::vector<std::pair<entity_key_t,body>>&,
  std::array<point_t,2>&,
  double);

void 
mpi_sort_unbalanced(
  std::vector<std::pair<entity_key_t,body>>&,
  int);

void 
mpi_branches_exchange_useful_positions(
    tree_topology_t& tree,
    std::vector<std::pair<entity_key_t,body>>&,
    std::vector<std::pair<point_t,point_t>>&,
    double);

void 
mpi_tree_traversal_graphviz(
  tree_topology_t&,
  std::array<point_t,2>&);

void 
mpi_compute_ghosts(
  tree_topology_t&,
  double smoothinglength,
  mpi_ghosts_t&,
  std::array<point_t,2>&);

void 
mpi_refresh_ghosts(
  tree_topology_t&,
  mpi_ghosts_t&,
  std::array<point_t,2>&);

void 
mpi_output_txt(
  std::vector<std::pair<entity_key_t,body>>&,
  int);

void 
computeAcceleration(
    point_t,
    point_t,
    double,
    point_t&,
    //point_t&,
    double*,
    double*);

bool 
MAC(
    branch_t*,
    branch_t*,
    double);

void 
tree_traversal_c2c(
    tree_topology_t&,
    branch_t*,
    branch_t*,
    point_t&,
    //point_t&,
    double*,
    double*,
    double&);

void 
sink_traversal_c2p(
  tree_topology_t&,
  branch_t*,
  point_t&,
  //point_t&,
  point_t&, 
  double*,
  double*,
  std::vector<body*>&,
  int&);

void 
tree_traversal_grav(
    tree_topology_t&,
    branch_t*,
    double&, 
    double&,
    int&);

void
tree_traversal_com(
    tree_topology_t&);



#if 0
void 
mpi_sort(
  std::vector<std::pair<entity_key_t,body>>&,
  std::vector<int>);
#endif

#if 0
void 
mpi_branches_exchange(
  tree_topology_t&);
#endif

#if 0
void 
mpi_branches_exchange_useful(
  tree_topology_t&,
  std::vector<std::pair<entity_key_t,body>>&,
  std::array<point_t,2>&,
  std::vector<std::pair<entity_key_t,entity_key_t>>&);
#endif 

#if 0
void 
mpi_gather_com(
  tree_topology_t&, 
  std::array<point_t,2>&,
  std::vector<std::pair<entity_key_t,entity_key_t>>&,
  std::vector<body_holder>&);
#endif 

#if 0
void 
mpi_gather_com_positions(
  tree_topology_t&, 
  std::array<point_t,2>&,
  std::vector<std::pair<point_t,point_t>>&,
  std::vector<body_holder>&);
#endif 

#if 0
void
traversal_COM(
  int,
  tree_topology_t&,
  std::array<point_t,2>&,
  branch_t *, 
  std::vector<body_holder>&,
  std::pair<point_t,point_t>&,
  int&);
#endif 

#if 0
bool
MAC(
    body_holder*,
    branch_t *,
    double,
    double);
#endif 

#if 0
void
traversal_COM_MAC_seq(
  tree_topology_t&,
  body_holder* bi,
  branch_t *, 
  std::vector<body_holder>&,
  double&,
  double&);
#endif



#endif // mpi_partition_h

