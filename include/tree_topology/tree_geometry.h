/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_topology_tree_geometry_h
#define flecsi_topology_tree_geometry_h

/*!
  \file tree_geometry.h
  \authors jloiseau@lanl.gov
  \date Initial file creation: Oct 9, 2018
 */

#include <map>
#include <unordered_map>
#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <functional>
#include <mutex>
#include <stack>
#include <math.h>
#include <float.h>

#include "flecsi/geometry/point.h"
#include "flecsi/concurrency/thread_pool.h"
#include "flecsi/data/storage.h"
#include "flecsi/data/data_client.h"
#include "flecsi/topology/index_space.h"


/*
#define np(X)                                                            \
 std::cout << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ \
           << ": " << #X << " = " << (X) << std::endl

#define hp(X)                                                            \
 std::cout << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ \
           << ": " << #X << " = " << std::hex << (X) << std::endl
*/

namespace flecsi {
namespace topology {

template<
  typename T,
  size_t D
>
struct tree_geometry{};

/*!
  \brief 1d geometry class for computing intersections and distances.
 */
template<
  typename T
>
struct tree_geometry<T, 1>
{
  using point_t = point__<T, 1>;
  using element_t = T;

  /*!
    Return true if point origin lies within the spheroid centered at center
    with radius.
   */
  static
  bool
  within(
    const point_t& origin,
    const point_t& center,
    element_t r1,
    element_t r2 = 0)
  {
    return distance(origin, center) <= r1;
  }

  /*!
    Return true if dist^2 < radius^2
   */
  static
  bool
  within_square(
    const point_t& origin,
    const point_t& center,
    element_t r1,
    element_t r2)
  {
    element_t dist = distance(origin,center);
    return dist*dist < (r1+r2)*(r1+r2);
  }

  /*!
    Return true if point origin lies within the box specified by min/max point.
   */
  static
  bool
  within_box(
    const point_t& min, 
    const point_t& max, 
    const point_t& origin,
    const element_t& r)
  {
    return origin[0] <= max[0] && origin[0] >= min[0];
  }

  static 
  bool
  intersects_box_box(
    const point_t& min_b1,
    const point_t& max_b1,
    const point_t& min_b2,
    const point_t& max_b2)
  {
    return 
      min_b1[0] <= max_b2[0] && min_b1[0] >= min_b2[0] || // b2 b1 b2 b1 
      max_b1[0] <= max_b2[0] && max_b1[0] >= min_b2[0] || // b1 b2 b1 b2 
      max_b1[0] >= max_b2[0] && min_b1[0] <= min_b2[0] || // b1 b2 b2 b1
      min_b1[0] >= min_b2[0] && max_b1[0] <= max_b2[0];   // b2 b1 b1 b2
  }

  // initial attempt to get this working, needs to be optimized

  static
  bool
  intersects(
    const point_t& origin,
    element_t size,
    const point_t& scale,
    const point_t& center,
    element_t radius)
  {
    return center[0] >= origin[0] - radius &&
           center[0] <= origin[0] + size * scale[0] + radius;
  }

  // Intersection of two spheres
  static 
  bool 
  intersects_sphere_sphere(
    const point_t& c1, 
    const element_t r1, 
    const point_t& c2, 
    const element_t r2)
  {
    return distance(c1,c2) <= r1+r2; 
  }

  // Intersection of sphere and box
  static 
  bool 
  intersects_sphere_box(
    const point_t& min, 
    const point_t& max, 
    const point_t& c, 
    const element_t r)
  {
    point_t x = point_t(std::max(min[0],std::min(c[0],max[0])));  
    element_t dist = distance(x,c); 
    return dist <= r;
  }

  static
  bool
  intersect_true(
    const point_t& origin,
    element_t size,
    const point_t& scale,
    const point_t& center,
    element_t radius)
  {
    return true;
  }

  static
  bool 
  within_mac(
    const point_t& p1, 
    const point_t& p2,
    const element_t radius,
    const element_t MAC)
  {
    return radius/distance(p1,p2) < MAC;
  }

  bool
  box_MAC(
    const point_t& position_source,
    const point_t& position_sink,
    const point_t& box_source_min,
    const point_t& box_source_max,
    double macangle)
  {
    double dmax = flecsi::distance(box_source_min,box_source_max);
    double disttoc = flecsi::distance(
        position_sink,position_source);
    return dmax/disttoc < macangle;
  }

};


/*!
  \brief 2d geometry class for computing intersections and distances.
 */
template<
  typename T
>
struct tree_geometry<T, 2>
{
  using point_t = point__<T, 2>;
  using element_t = T;

  /*!
    Return true if point origin lies within the spheroid centered at center
    with radius.
   */
  static
  bool
  within(
    const point_t& origin,
    const point_t& center,
    element_t r1,
    element_t r2 = 0)
  {
    return distance(origin, center) <= r1;
  }

  /*!
    Return true if dist^2 < radius^2
   */
  static
  bool
  within_square(
    const point_t& origin,
    const point_t& center,
    element_t r1,
    element_t r2)
  {
    element_t dist = distance(origin,center);
    return dist*dist < (r1+r2)*(r1+r2);
  }

  /*!
    Return true if point origin lies within the box specified by min/max point.
   */
  static
  bool
  within_box(
    const point_t& min,
    const point_t& max,
    const point_t& origin,
    const element_t& r)
  {
    return origin[0] <= max[0] && origin[0] >= min[0] &&
           origin[1] <= max[1] && origin[1] >= min[1];
  }

  // initial attempt to get this working, needs to be optimized

  /*!
    Spheroid/box intersection test.
   */
  static
  bool
  intersects(
    const point_t& origin,
    element_t size,
    const point_t& scale,
    const point_t& center,
    element_t radius)
  {
    return center[0] >= origin[0] - radius &&
           center[0] <= origin[0] + size * scale[0] + radius &&
           center[1] >= origin[1] - radius &&
           center[1] <= origin[1] + size * scale[1] + radius;
  }

  static
  bool 
  within_mac(
    const point_t& p1, 
    const point_t& p2,
    const element_t radius,
    const element_t MAC)
  {
    return 2*asin(radius/distance(p1,p2)) < MAC;
  }

  static 
  bool
  intersects_box_box(
    const point_t& min_b1,
    const point_t& max_b1,
    const point_t& min_b2,
    const point_t& max_b2)
  {
    return 
      (min_b1[0] <= max_b2[0] && min_b1[0] >= min_b2[0] || // b2 b1 b2 b1 
       max_b1[0] <= max_b2[0] && max_b1[0] >= min_b2[0] || // b1 b2 b1 b2 
       max_b1[0] >= max_b2[0] && min_b1[0] <= min_b2[0] || // b1 b2 b2 b1
       min_b1[0] >= min_b2[0] && max_b1[0] <= max_b2[0])   // b2 b1 b1 b2
      &&
      (min_b1[1] <= max_b2[1] && min_b1[1] >= min_b2[1] || // b2 b1 b2 b1 
       max_b1[1] <= max_b2[1] && max_b1[1] >= min_b2[1] || // b1 b2 b1 b2 
       max_b1[1] >= max_b2[1] && min_b1[1] <= min_b2[1] || // b1 b2 b2 b1
       min_b1[1] >= min_b2[1] && max_b1[1] <= max_b2[1]);  // b2 b1 b1 b2
  }

  // Intersection of two spheres
  static 
  bool 
  intersects_sphere_sphere(
    const point_t& c1, 
    const element_t r1, 
    const point_t& c2, 
    const element_t r2)
  {
    return distance(c1,c2) <= r1+r2; 
  }

  // Intersection of sphere and box
  static 
  bool 
  intersects_sphere_box(
    const point_t& min, 
    const point_t& max, 
    const point_t& c, 
    const element_t r)
  {
    point_t x = point_t(
        std::max(min[0],std::min(c[0],max[0])),
        std::max(min[1],std::min(c[1],max[1])));  
    element_t dist = distance(x,c); 
    return dist <= r;
  }


  static
  bool
  intersect_true(
    const point_t& origin,
    element_t size,
    const point_t& scale,
    const point_t& center,
    element_t radius)
  {
    return true;
  }

  static
  bool
  box_MAC(
    const point_t& position_source,
    const point_t& position_sink,
    const point_t& box_source_min,
    const point_t& box_source_max,
    double macangle)
  {
    double dmax = flecsi::distance(box_source_min,box_source_max);
    double disttoc = flecsi::distance(
        position_sink,position_source);
    return dmax/disttoc < macangle;
  }

};


/*!
  \brief 1d geometry class for computing intersections and distances.
 */
template<
  typename T
>
struct tree_geometry<T, 3>
{
  using point_t = point__<T, 3>;
  using element_t = T;

  /*!
    Return true if point origin lies within the spheroid centered at center
    with radius.
   */
  static
  bool
  within(
    const point_t& origin,
    const point_t& center,
    element_t r1,
    element_t r2 = 0)
  {
    return distance(origin, center) <= r1;
  }

  /*!
    Return true if dist^2 < radius^2
   */
  static
  bool
  within_square(
    const point_t& origin,
    const point_t& center,
    element_t r1,
    element_t r2)
  {
    element_t dist = distance(origin,center);
    return dist*dist < (r1+r2)*(r1+r2);
  }
  
  static
  bool 
  within_mac(
    const point_t& p1, 
    const point_t& p2,
    const element_t radius,
    const element_t MAC)
  {
    return 2*asin(radius/distance(p1,p2)) < MAC;
  }

  static 
  bool
  intersects_box_box(
    const point_t& min_b1,
    const point_t& max_b1,
    const point_t& min_b2,
    const point_t& max_b2)
  {
    return 
      (min_b1[0] <= max_b2[0] && min_b1[0] >= min_b2[0] || // b2 b1 b2 b1 
       max_b1[0] <= max_b2[0] && max_b1[0] >= min_b2[0] || // b1 b2 b1 b2 
       max_b1[0] >= max_b2[0] && min_b1[0] <= min_b2[0] || // b1 b2 b2 b1
       min_b1[0] >= min_b2[0] && max_b1[0] <= max_b2[0])   // b2 b1 b1 b2
      &&
      (min_b1[1] <= max_b2[1] && min_b1[1] >= min_b2[1] || // b2 b1 b2 b1 
       max_b1[1] <= max_b2[1] && max_b1[1] >= min_b2[1] || // b1 b2 b1 b2 
       max_b1[1] >= max_b2[1] && min_b1[1] <= min_b2[1] || // b1 b2 b2 b1
       min_b1[1] >= min_b2[1] && max_b1[1] <= max_b2[1])   // b2 b1 b1 b2
      &&
      (min_b1[2] <= max_b2[2] && min_b1[2] >= min_b2[2] || // b2 b1 b2 b1 
       max_b1[2] <= max_b2[2] && max_b1[2] >= min_b2[2] || // b1 b2 b1 b2 
       max_b1[2] >= max_b2[2] && min_b1[2] <= min_b2[2] || // b1 b2 b2 b1
       min_b1[2] >= min_b2[2] && max_b1[2] <= max_b2[2]);  // b2 b1 b1 b2
  }


  /*!
    Return true if point origin lies within the box specified by min/max point.
   */
  static
  bool
  within_box(
    const point_t& min,
    const point_t& max,
    const point_t& origin,
    const element_t& r)
  {
    return origin[0] <= max[0] && origin[0] > min[0] &&
           origin[1] <= max[1] && origin[1] > min[1] &&
           origin[2] <= max[2] && origin[2] > min[2];
  }



  /*!
    Spheroid/box intersection test.
   */
  static
  bool
  intersects(
    const point_t& origin,
    element_t size,
    const point_t& scale,
    const point_t& center,
    element_t radius)
  {
    return center[0] >= origin[0] - radius &&
           center[0] <= origin[0] + size * scale[0] + radius &&
           center[1] >= origin[1] - radius &&
           center[1] <= origin[1] + size * scale[1] + radius &&
           center[2] >= origin[2] - radius &&
           center[2] <= origin[2] + size * scale[2] + radius;
  }
  
  // Intersection of two spheres
  static 
  bool 
  intersects_sphere_sphere(
    const point_t& c1, 
    const element_t r1, 
    const point_t& c2, 
    const element_t r2)
  {
    return (c2[0]-c1[0])*(c2[0]-c1[0])+
      (c2[1]-c1[1])*(c2[1]-c1[1])+
      (c2[2]-c1[2])*(c2[2]-c1[2])
        < (r1+r2)*(r1+r2);
  }



  // Intersection of sphere and box
  static 
  bool 
  intersects_sphere_box(
    const point_t& min, 
    const point_t& max, 
    const point_t& c, 
    const element_t r)
  {
    point_t x = point_t(
        std::max(min[0],std::min(c[0],max[0])),
        std::max(min[1],std::min(c[1],max[1])), 
        std::max(min[2],std::min(c[2],max[2]))); 
    element_t dist = distance(x,c); 
    return dist <= r;
  }

  static
  bool
  intersect_true(
    const point_t& origin,
    element_t size,
    const point_t& scale,
    const point_t& center,
    element_t radius)
  {
    return true;
  }

  static
  bool
  box_MAC(
    const point_t& position_source,
    const point_t& position_sink,
    const point_t& box_source_min,
    const point_t& box_source_max,
    double macangle)
  {
    double dmax = flecsi::distance(box_source_min,box_source_max);
    double disttoc = flecsi::distance(
        position_sink,position_source);
    return dmax/disttoc < macangle;
  }

};
} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_tree_geometry_h
