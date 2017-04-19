#ifndef io_h
#define io_h 

#include <cstdlib>
#include <iostream>
#include <vector>

#include <hdf5.h>

#include "tree.h"

namespace io{
  void inputDataTxt(std::vector<body*>&, const char *, tree_topology_t&); 
  void inputDataHDF5(std::vector<body*>&, const char *, tree_topology_t&); 
  void inputDataTxtRange(std::vector<body*>&,int&,int&,int,int,const char *);

  void outputDataTxt(); 
  void outputDataHDF5(std::vector<body*>&,const char*, int, double); 
} // namespace io

#endif // io_h
