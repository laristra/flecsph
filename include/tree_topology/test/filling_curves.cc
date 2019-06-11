#include <cinchdevel.h>
#include <cinchtest.h>

#include <cmath>
#include <iostream>
#include <mpi.h>

#include "../filling_curve.h"

using namespace std;
using namespace flecsi;

namespace flecsi {
namespace execution {
void driver(int argc, char *argv[]) {}
} // namespace execution
} // namespace flecsi

template <typename T, size_t D>
bool operator==(const point__<T, D> &a, const point__<T, D> &b) {
  for (size_t d = 0; d < D; ++d) {
    if (a[d] != b[d])
      return false;
  }
  return true;
}

using point_t = point__<double, 3>;
using range_t = std::array<point_t, 2>;
using hc = hilbert_curve_u<3, uint64_t>;
using mc = morton_curve_u<3, uint64_t>;

using point_2d = point__<double,2>;
using range_2d = std::array<point_2d,2>;
using hc_2d = hilbert_curve_u<2,uint64_t>;
using mc_2d = morton_curve_u<2,uint64_t>;

TEST(hilbert,sanity) {

  using namespace flecsi;

  range_t range;
  range[0] = {0, 0, 0};
  range[1] = {1, 1, 1};
  point_t p1 = {0.25, 0.25, 0.25};

  std::cout << "Hilbert TEST " << hc::max_depth() << std::endl;

  hc::set_range(range);
  hc hc1;
  hc hc2(p1);
  hc hc3 = hc::min();
  hc hc4 = hc::max();
  hc hc5 = hc::root();
  std::cout << "Default: " << hc1 << std::endl;
  std::cout << "Pt:rg  : " << hc2 << std::endl;
  std::cout << "Min    : " << hc3 << std::endl;
  std::cout << "Max    : " << hc4 << std::endl;
  std::cout << "root   : " << hc5 << std::endl;
  ASSERT_TRUE(1 == hc5);

  while(hc4 != hc5) {
    hc4.pop();
  }
  ASSERT_TRUE(hc5 == hc4);

}

#if 0 
TEST(hilbert,rnd_2d){
  using namespace flecsi;
  // Test the generation 2D
  range_2d rge;
  rge[0] = {0,0};
  rge[1] = {1,1};
  hc_2d::set_range(rge);
  std::array<point_2d,4> pts = {
    point_2d{.25,.25},point_2d{.25,.5},point_2d{.5,.5},point_2d{.5,.25}};
  std::array<hc_2d,4> hcs_2d;

  for(int i = 0 ; i < 4; ++i){
    hcs_2d[i] = hc_2d(pts[i]);
    point_2d inv;
    hcs_2d[i].coordinates(inv);
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
  hc::set_range(range);
  std::array<point_t,8> points = {
    point_t{.25,.25,.25},point_t{.25,.25,.5},point_t{.25,.5,.5},point_t{.25,.5,.25},
    point_t{.5,.5,.25},point_t{.5,.5,.5},point_t{.5,.25,.5},point_t{.5,.25,.25}};
  std::array<hc,8> hcs;

  for(int i = 0 ; i < 8; ++i){
    hcs[i] = hc(points[i]);
    point_t inv;
    hcs[i].coordinates(inv);
    double dist = distance(points[i],inv);
    std::cout << points[i] <<" "<< hcs[i] << " = "<<inv<<std::endl;
    ASSERT_TRUE(dist<1.0e-4);
  }

  // rnd
  for(int i = 0 ; i < 20; ++i){
    point_t pt(
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX);
    point_t inv;
    hc h(pt);
    h.coordinates(inv);
    double dist = distance(pt,inv);
    std::cout << pt <<" = "<< h << " = "<<inv<<std::endl;
    ASSERT_TRUE(dist<1.0e-4);
  }
} // TEST
#endif 

TEST(morton,sanity) {
  range_t range;
  range[0] = {-1, -1, -1};
  range[1] = {1, 1, 1};
  point_t p1 = {0, 0, 0};

  std::cout <<" Morton TEST "<< hc::max_depth() << std::endl;

  mc::set_range(range);
  mc hc1;
  mc hc2(p1);
  mc hc3 = mc::min();
  mc hc4 = mc::max();
  mc hc5 = mc::root();
  std::cout << "Default: " << hc1 << std::endl;
  std::cout << "Pt:rg  : " << hc2 << std::endl;
  std::cout << "Min    : " << hc3 << std::endl;
  std::cout << "Max    : " << hc4 << std::endl;
  std::cout << "root   : " << hc5 << std::endl;
  ASSERT_TRUE(1 == hc5);

  while(hc4 != hc5) {
    hc4.pop();
  }
  ASSERT_TRUE(hc5 == hc4);
}

TEST(morton,rnd_2d){
  using namespace flecsi;
  // Test the generation 2d
  range_2d rge;
  rge[0] = {0,0};
  rge[1] = {1,1};
  mc_2d::set_range(rge);
  std::array<point_2d,4> pts = {
    point_2d{.25,.25},point_2d{.5,.25},point_2d{.25,.5},point_2d{.5,.5}};
  std::array<mc_2d,4> mcs_2d;

  for(int i = 0 ; i < 4; ++i){
    mcs_2d[i] = mc_2d(pts[i]);
    point_2d inv;
    mcs_2d[i].coordinates(inv);
    double dist = distance(pts[i],inv);
    std::cout << pts[i] <<" "<< mcs_2d[i] << " = "<<inv<<std::endl;
    ASSERT_TRUE(dist<1.0e-4);
  }

  // rnd
  for(int i = 0 ; i < 20; ++i){
    point_2d pt(
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX);
    point_2d inv;
    mc_2d h(pt);
    h.coordinates(inv);
    double dist = distance(pt,inv);
    std::cout << pt <<" = "<< h << " = "<<inv<<std::endl;
    ASSERT_TRUE(dist<1.0e-4);
  }
}

TEST(morton,rnd_3d){
  using namespace flecsi;
  range_t range;

  // Test the generation
  range[0] = {0,0,0};
  range[1] = {1,1,1};
  mc::set_range(range);
  std::array<point_t,8> points = {
    point_t{.25,.25,.25},point_t{.5,.25,.25},point_t{.25,.5,.25},point_t{.5,.5,.25},
    point_t{.25,.25,.5},point_t{.5,.25,.5},point_t{.25,.5,.5},point_t{.5,.5,.5}};
  std::array<mc,8> mcs;

  for(int i = 0 ; i < 8; ++i){
    mcs[i] = mc(points[i]);
    point_t inv;
    mcs[i].coordinates(inv);
    double dist = distance(points[i],inv);
    std::cout << points[i] <<" "<< mcs[i] << " = "<<inv<<std::endl;
    ASSERT_TRUE(dist<1.0e-4);
  }

  // rnd
  for(int i = 0 ; i < 20; ++i){
    point_t pt(
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX,
      (double)rand()/(double)RAND_MAX);
    point_t inv;
    mc h(pt);
    h.coordinates(inv);
    double dist = distance(pt,inv);
    std::cout << pt <<" = "<< h << " = "<<inv<<std::endl;
    ASSERT_TRUE(dist<1.0e-4);
  }
}
