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

void mpi_compute_range(
    std::vector<std::pair<entity_key_t,body>>&,
    std::array<point_t,2>&);

void mpi_sort(
    std::vector<std::pair<entity_key_t,body>>&,
    std::vector<int>);

void mpi_branches_exchange(
  tree_topology_t&);

void mpi_tree_traversal_graphviz(
    tree_topology_t&,
    std::array<point_t,2>&);

void mpi_gather_ghosts(
    tree_topology_t&,
    double smoothinglength,
    std::vector<body>&);

void mpi_output_txt(
    std::vector<std::pair<entity_key_t,body>>&,
    std::vector<int>&,
    int
    );

#endif // mpi_partition_h
