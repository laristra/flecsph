#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "kernel.h"

using namespace std;
using namespace flecsi;
using namespace topology;
using namespace kernel;

namespace flecsi{
namespace execution{
  void driver(int argc, char* argv[]){
  }
}
}


TEST(kernel, cubic_spline) {

  // Generate particles, write to file, read and compare 
  double h = .05;
  double step = h/100.;
  double max = h * 6.;
  size_t n = max / step;
  
  double current = -3.*h;

  point_t p;

  char name[] = "output_kernel.csv";

  // Output 
  FILE * output = fopen(name,"w");

  for(size_t i = 0; i < n; ++i){
    p[0] = current;
    fprintf(output,"%.4f;%.4f;%.4f\n",
      current,
      cubic_spline_kernel(fabs(current),h),
      cubic_spline_gradKernel(p,h)[0]);
    current += step; 
  }

  fclose(output);
}