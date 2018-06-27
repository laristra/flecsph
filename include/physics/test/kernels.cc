#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "kernels.h"

using namespace std;
using namespace flecsi;
using namespace topology;
using namespace kernels;

namespace flecsi{
namespace execution{
  void driver(int argc, char* argv[]){
  }
}
}


TEST(kernel, cubic_spline) {

  // Generate particles, write to file, read and compare 
  double h = 1.;
  double step = 0.01;
  double max = h * 6.;
  size_t n = max / step;
  
  double current = -3.*h;

  point_t p;

  char name[] = "cubic_spline.csv";

  // Output 
  FILE * output = fopen(name,"w");

  for(size_t i = 0; i < n; ++i){
    p[0] = current;
    fprintf(output,"%.4f;%.4f;%.4f\n",
      current,
      cubic_spline(fabs(current),h),
      gradient_cubic_spline(p,h)[0]);
    current += step; 
  }

  fclose(output);
}

TEST(kernel, gaussian) {

  // Generate particles, write to file, read and compare 
  double h = 1.;
  double step = 0.01;
  double max = h * 6.;
  size_t n = max / step;
  
  double current = -3.*h;

  point_t p;

  char name[] = "gaussian.csv";

  // Output 
  FILE * output = fopen(name,"w");

  for(size_t i = 0; i < n; ++i){
    p[0] = current;
    fprintf(output,"%.4f;%.4f;%.4f\n",
      current,
      gaussian(fabs(current),h),
      gradient_gaussian(p,h)[0]);
    current += step; 
  }

  fclose(output);
}

TEST(kernel, quintic_spline) {

  // Generate particles, write to file, read and compare 
  double h = 1.;
  double step = 0.01;
  double max = h * 6.;
  size_t n = max / step;
  
  double current = -3.*h;

  point_t p;

  char name[] = "quintic_spline.csv";

  // Output 
  FILE * output = fopen(name,"w");

  for(size_t i = 0; i < n; ++i){
    p[0] = current;
    fprintf(output,"%.4f;%.4f;%.4f\n",
      current,
      quintic_spline(fabs(current),h),
      gradient_quintic_spline(p,h)[0]);
    current += step; 
  }

  fclose(output);
}

TEST(kernel, wendland_quintic) {

  // Generate particles, write to file, read and compare 
  double h = 1.;
  double step = 0.01;
  double max = h * 6.;
  size_t n = max / step;
  
  double current = -3.*h;

  point_t p;

  char name[] = "wendland_quintic.csv";

  // Output 
  FILE * output = fopen(name,"w");

  for(size_t i = 0; i < n; ++i){
    p[0] = current;
    fprintf(output,"%.4f;%.4f;%.4f\n",
      current,
      wendland_quintic(fabs(current),h),
      gradient_wendland_quintic(p,h)[0]);
    current += step; 
  }

  fclose(output);
}