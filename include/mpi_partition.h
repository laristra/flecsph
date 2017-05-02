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

void mpi_compute_range(
    std::vector<std::pair<entity_key_t,body>>&,
    std::array<point_t,2>&);

void mpi_sort_unbalanced(
    std::vector<std::pair<entity_key_t,body>>&,
    int);

void mpi_sort(
    std::vector<std::pair<entity_key_t,body>>&,
    std::vector<int>);

void mpi_branches_exchange(
  tree_topology_t&);

void mpi_branches_exchange_useful(
    tree_topology_t&,
    std::vector<std::pair<entity_key_t,body>>&,
    std::array<point_t,2>&);

void mpi_tree_traversal_graphviz(
    tree_topology_t&,
    std::array<point_t,2>&);

void mpi_compute_ghosts(
    tree_topology_t&,
    double smoothinglength,
    mpi_ghosts_t&);

void mpi_refresh_ghosts(
    tree_topology_t&,
    mpi_ghosts_t&);

void mpi_output_txt(
    std::vector<std::pair<entity_key_t,body>>&,
    int
    );

#endif // mpi_partition_h
