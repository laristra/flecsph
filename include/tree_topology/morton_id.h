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

#ifndef flecsi_topology_morton_id_h
#define flecsi_topology_morton_id_h

/*!
  \file morton_id.h
  \authors jloiseau@lanl.gov
  \date October 9, 2018
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
#include <math.h>
#include <float.h>

#include "flecsi/geometry/point.h"
#include "flecsi/concurrency/thread_pool.h"
#include "flecsi/data/storage.h"
#include "flecsi/data/data_client.h"
#include "flecsi/topology/index_space.h"
 

namespace flecsi{
namespace topology{

/*!
  This class implements a hashed/Morton-style id that can be
  parameterized on arbitrary dimension D and integer type T.
 */
template<
  typename T,
  size_t D
>
class morton_id
{
public:

  using int_t = T;
  
  static const size_t dimension = D;
  static constexpr size_t bits = sizeof(int_t) * 8;
  static constexpr size_t max_depth = (bits - 1)/dimension;

  morton_id()
  : id_(0)
  {}

  /*!
    Construct the id from an array of dimensions and range for each
    dimension. The specified depth may be less than the max allowed depth for
    the id.
   */
  template<
    typename S
  >
  morton_id(
    const std::array<point__<S, dimension>, 2>& range,
    const point__<S, dimension>& p,
    size_t depth)
  : id_(int_t(1) << depth * dimension + (bits - 1) % dimension)
  {
    std::array<int_t, dimension> coords;

    for(size_t i = 0; i < dimension; ++i)
    {
      S min = range[0][i];
      S scale = range[1][i] - min;
      coords[i] = (p[i] - min)/scale * (int_t(1) << (bits - 1)/dimension);
    }

    size_t k = 0;
    for(size_t i = max_depth - depth; i < max_depth; ++i)
    {
      for(size_t j = 0; j < dimension; ++j)
      {
        int_t bit = (coords[j] & int_t(1) << i) >> i;
        id_ |= bit << (k * dimension + j);
      }
      ++k;
    }
  }

  /*!
    Construct the id from an array of dimensions and range for each
    dimension. Construct with the max depth available.
   */
  template<
    typename S
  >
  morton_id(
    const std::array<point__<S, dimension>, 2>& range,
    const point__<S, dimension>& p)
  : morton_id(range,p,max_depth)
  {}


  constexpr morton_id(const morton_id& bid)
  : id_(bid.id_)
  {}

  /*!
    Get the root id (depth 0).
   */
  static
  constexpr
  morton_id
  root()
  {
    return morton_id(int_t(1) << (bits - 1) % dimension);
  }

  /*!
    Get the null id.
   */
  static
  constexpr
  morton_id
  null()
  {
    return morton_id(0);
  }

  /*!
    Check if id is null.
   */
  constexpr
  bool
  is_null() const
  {
    return id_ == int_t(0);
  }

  /*!
    Find the depth of this id.
   */
  size_t
  depth() const
  {
    int_t id = id_;
    size_t d = 0;

    while(id >>= dimension)
    {
      ++d;
    }

    return d;
  }

  morton_id&
  operator=(
    const morton_id& bid
  )
  {
    id_ = bid.id_;
    return *this;
  }

  constexpr
  bool
  operator==(
    const morton_id& bid
  ) const
  {
    return id_ == bid.id_;
  }

  constexpr
  bool
  operator!=(
    const morton_id& bid
  ) const
  {
    return id_ != bid.id_;
  }

  /*!
    Push bits onto the end of this id.
   */
  void push(int_t bits)
  {
    assert(bits < int_t(1) << dimension);

    id_ <<= dimension;
    id_ |= bits;
  }

  /*!
    Pop the bits of greatest depth off this id.
   */
  void pop()
  {
    assert(depth() > 0);
    id_ >>= dimension;
  }

  /*!
    Pop the depth d bits from the end of this this id.
   */
  void pop(
    size_t d
  )
  {
    assert(d >= depth());
    id_ >>= d * dimension;
  }

  /*!
    Return the parent of this id (depth - 1)
   */
  constexpr
  morton_id
  parent() const
  {
    return morton_id(id_ >> dimension);
  }

  /*!
    Truncate (repeatedly pop) this id until it of depth to_depth.
   */
  void
  truncate(
    size_t to_depth
  )
  {
    size_t d = depth();

    if(d < to_depth)
    {
      return;
    }

    id_ >>= (d - to_depth) * dimension;
  }

  void
  output_(
    std::ostream& ostr
  ) const
  {
    constexpr int_t mask = ((int_t(1) << dimension) - 1) << bits - dimension;

    size_t d = max_depth;

    int_t id = id_;

    while((id & mask) == int_t(0))
    {
      --d;
      id <<= dimension;
    }

    if(d == 0)
    {
      ostr << "1";
      return;
    }
    ostr << "1";
    //ostr<<std::bitset<64>(id);
    id <<= 1 + (bits - 1) % dimension;

    for(size_t i = 1; i <= d; ++i)
    {
      int_t val = (id & mask) >> (bits - dimension);
      ostr << std::bitset<D>(val);
      id <<= dimension;
    }
  }

  int_t
  value_() const
  {
    return id_;
  }

  bool
  operator<(
    const morton_id& bid
  ) const
  {
    return id_ < bid.id_;
  }

  /*!
    Convert this id to coordinates in range.
   */
  template<
    typename S
  >
  void
  coordinates(
    const std::array<point__<S, dimension>, 2>& range,
    point__<S, dimension>& p) const
  {
    std::array<int_t, dimension> coords;
    coords.fill(int_t(0));

    int_t id = id_;
    size_t d = 0;

    while(id >> dimension != int_t(0))
    {
      for(size_t j = 0; j < dimension; ++j)
      {
        coords[j] |= (((int_t(1) << j) & id) >> j) << d;
      }

      id >>= dimension;
      ++d;
    }

    constexpr int_t m = (int_t(1) << max_depth) - 1;

    for(size_t j = 0; j < dimension; ++j)
    {
      S min = range[0][j];
      S scale = range[1][j] - min;

      coords[j] <<= max_depth - d;
      p[j] = min + scale * S(coords[j])/m;
    }
  }

  /**
   * @brief Compute the range of a branch from its key 
   * The space is recursively decomposed regarding the dimension 
   */
  template<
    typename S
  >
  std::array<point__<S,dimension>, 2>
  range(
      const std::array<point__<S,dimension>, 2>& range)
  { 
    // The result range 
    std::array<point__<S,dimension>, 2> result;
    result[0] = range[0]; 
    result[1] = range[1]; 
    // Copy the key 
    int_t tmp = id_; 
    int_t root = morton_id::root().id_; 
    
    // Extract x,y and z 
    std::array<int_t, dimension> coords;
    coords.fill(int_t(0));

    int_t id = id_;
    size_t d = 0;

    while(id != root)
    {
      for(size_t j = 0; j < dimension; ++j)
      {
        coords[j] |= (((int_t(1) << j) & id) >> j) << d;
      }

      id >>= dimension;
      ++d;
    }

    std::cout<<"depth="<<d<<std::endl;
   
    for(size_t i = 0 ; i < dimension ; ++i)
    {
      // apply the reduction 
      for(size_t j = d ; j > 0; --j)
      {
        double nu = (result[0][i]+result[1][i])/2.;
        if(coords[i] & (int_t(1)<<j-1))
        {
          result[0][i] = nu; 
        }else{
          result[1][i] = nu;
        }
      }
    }
    return result; 
  }

private:
  int_t id_;

  constexpr
  morton_id(
    int_t id
  )
  : id_(id)
  {}
};

// output for morton id 
template<
  typename T, 
  size_t D> 
std::ostream&
operator<<(
  std::ostream& ostr,
  const morton_id<T,D>& k)
{
  k.output_(ostr);
  return ostr;
}

} // namespace topology 
} // namespace flecsi

#endif // flecsi_topology_morton_id_h
