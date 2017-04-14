#ifndef io_h
#define io_h 

#include <cstdlib>
#include <iostream>
#include <vector>

#include "tree.h"

namespace io{
  void inputDataTxt(std::vector<body*>&, const char *, tree_topology_t&); 
} // namespace io

#endif // io_h
