#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <log.h>
#include <mpi.h>

#include "../filling_curve.h"

using namespace std;
using namespace flecsi;

namespace flecsi {
namespace execution {
void
driver(int, char **) {}
} // namespace execution
} // namespace flecsi

template<typename T, size_t D>
bool
operator==(const space_vector_u<T, D> & a, const space_vector_u<T, D> & b) {
  for(size_t d = 0; d < D; ++d) {
    if(a[d] != b[d])
      return false;
  }
  return true;
}

using point_t = space_vector_u<double, 3>;
using range_t = std::array<point_t, 2>;
using hc = hilbert_curve_u<3, uint64_t>;
using mc = morton_curve_u<3, uint64_t>;

using point_2d = space_vector_u<double, 2>;
using range_2d = std::array<point_2d, 2>;
using hc_2d = hilbert_curve_u<2, uint64_t>;
using mc_2d = morton_curve_u<2, uint64_t>;

TEST(hilbert, sanity) {

  using namespace flecsi;

  range_t range;
  range[0] = {0, 0, 0};
  range[1] = {1, 1, 1};
  point_t p1 = {0.25, 0.25, 0.25};

  std::cout << "Hilbert TEST " << hc::max_depth() << std::endl;

  hc hc1;
  hc hc2(range, p1);
  hc hc3 = hc::min();
  hc hc4 = hc::max();
  hc hc5 = hc::root();
  std::cout << "Default: " << hc1 << std::endl;
  std::cout << "Pt:rg  : " << hc2 << std::endl;
  std::cout << "Min    : " << hc3 << std::endl;
  std::cout << "Max    : " << hc4 << std::endl;
  std::cout << "root   : " << hc5 << std::endl;
  ASSERT_TRUE(1 == hc5.value());

  while(hc4 != hc5) {
    hc4.pop();
  }
  ASSERT_TRUE(hc5 == hc4);
}

TEST(hilbert,rnd_2d){
  using namespace flecsi;
  // Test the generation 2D
  range_2d rge;
  rge[0] = {0,0};
  rge[1] = {1,1};
  const int npoints = 4; 
  std::array<point_2d,npoints> pts = {
    point_2d{0.,0.}, 
    point_2d{0.,1.},
    point_2d{1.,0.},
    point_2d{1.,1.},
  };
  std::array<hc_2d,4> hcs_2d;

  for(int i = 0 ; i < 4; ++i){
    hcs_2d[i] = hc_2d(rge,pts[i]);
    point_2d inv;
    hcs_2d[i].coordinates(rge,inv);
    double dist = distance(pts[i],inv);
    std::cout << pts[i] <<" "<< hcs_2d[i] << " = "<<inv<<std::endl;
    ASSERT_TRUE(dist<1.0e-4);
  }
}
  
TEST(hilbert,rnd_3d){
  using namespace flecsi;
  // Test the generation
  range_t range;

  range[0] = {0,0,0};
  range[1] = {1,1,1};
  const int npoints = 8; 
  std::array<point_t, npoints> points = {
    point_t{0.,0.,0.}, 
    point_t{0.,0.,1.},
    point_t{0.,1.,0.}, 
    point_t{0.,1.,1.}, 
    point_t{1.,0.,0.},
    point_t{1.,0.,1.}, 
    point_t{1.,1.,0.}, 
    point_t{1.,1.,1.}};
  std::array<hc,npoints> hcs;

  for(int i = 0 ; i < npoints; ++i){
    hcs[i] = hc(range,points[i]);
    point_t inv;
    hcs[i].coordinates(range,inv);
    //double dist = distance(points[i],inv);
    std::cout << points[i] <<" "<< hcs[i] << " = "<<inv<<std::endl;
    //ASSERT_TRUE(dist<1.0e-4);
  }

  // rnd
  for(int i = 0 ; i < 20; ++i){
    point_t pt(
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX);
    point_t inv;
    hc h(range,pt);
    h.coordinates(range,inv);
    //double dist = distance(pt,inv);
    std::cout << pt <<" = "<< h << " = "<<inv<<std::endl;
    //ASSERT_TRUE(dist<1.0e-4);
  }
} // TEST

TEST(morton, sanity) {
  range_t range;
  range[0] = {-1, -1, -1};
  range[1] = {1, 1, 1};
  point_t p1 = {0, 0, 0};

  std::cout << " Morton TEST " << hc::max_depth() << std::endl;

  // mc::set_range(range);
  mc hc1;
  mc hc2(range, p1);
  mc hc3 = mc::min();
  mc hc4 = mc::max();
  mc hc5 = mc::root();
  std::cout << "Default: " << hc1 << std::endl;
  std::cout << "Pt:rg  : " << hc2 << std::endl;
  std::cout << "Min    : " << hc3 << std::endl;
  std::cout << "Max    : " << hc4 << std::endl;
  std::cout << "root   : " << hc5 << std::endl;
  ASSERT_TRUE(1 == hc5.value());

  while(hc4 != hc5) {
    hc4.pop();
  }
  ASSERT_TRUE(hc5 == hc4);
}


TEST(morton, rnd_2d) {
  using namespace flecsi;
  // Test the generation 2d
  range_2d rge;
  rge[0] = {0., 0.};
  rge[1] = {1., 1.};
  const int npoints = 4; 
  std::array<point_2d, npoints> pts = {
    point_2d{0.,0.},
    point_2d{0.,1.},
    point_2d{1.,0.},
    point_2d{1.,1.}
  };
  std::array<mc_2d, npoints> mcs_2d;

  for(int i = 0; i < npoints; ++i) {
    mcs_2d[i] = mc_2d(rge, pts[i]);
    point_2d inv;
    mcs_2d[i].coordinates(rge, inv);
    double dist = distance(pts[i], inv);
    std::cout << pts[i] << " " << mcs_2d[i] << " = " << inv << std::endl;
    ASSERT_TRUE(dist < 1.0e-4);
  }

  // rnd
  for(int i = 0; i < 20; ++i) {
    point_2d pt(
      (double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX);
    point_2d inv;
    mc_2d h(rge, pt);
    h.coordinates(rge, inv);
    double dist = distance(pt, inv);
    std::cout << pt << " = " << h << " = " << inv << std::endl;
    ASSERT_TRUE(dist < 1.0e-4);
  }
}

TEST(morton, rnd_3d) {
  using namespace flecsi;
  range_t range;

  // Test the generation
  range[0] = {0, 0, 0};
  range[1] = {1, 1, 1};
  const int npoints = 8; 
  std::array<point_t, npoints> points = {
    point_t{0.,0.,0.}, 
    point_t{0.,0.,1.},
    point_t{0.,1.,0.}, 
    point_t{0.,1.,1.}, 
    point_t{1.,0.,0.},
    point_t{1.,0.,1.}, 
    point_t{1.,1.,0.}, 
    point_t{1.,1.,1.}};
  std::array<mc, 8> mcs;

  for(int i = 0; i < npoints; ++i) {
    mcs[i] = mc(range, points[i]);
    point_t inv;
    mcs[i].coordinates(range, inv);
    double dist = distance(points[i], inv);
    std::cout << points[i] << " " << mcs[i] << " = " << inv << std::endl;
    ASSERT_TRUE(dist < 1.0e-4);
  }

  // rnd
  for(int i = 0; i < 20; ++i) {
    point_t pt((double)rand() / (double)RAND_MAX,
      (double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX);
    point_t inv;
    mc h(range, pt);
    h.coordinates(range, inv);
    double dist = distance(pt, inv);
    std::cout << pt << " = " << h << " = " << inv << std::endl;
    ASSERT_TRUE(dist < 1.0e-4);
  }
}
