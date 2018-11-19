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

#ifndef flecsi_topology_hilbert_id_h
#define flecsi_topology_hilbert_id_h

/*!
  \file hilbert_id.h
  \authors jloiseau@lanl.gov
  \date November 16, 2018
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
class hilbert_id
{
public:

  using int_t = T;

  static const size_t dimension = D;
  static constexpr size_t bits = sizeof(int_t) * 8;
  static constexpr size_t max_depth = bits/dimension-1;

  hilbert_id()
  : id_(0)
  {}


  void rotation(int_t& n,
      std::array<int_t,dimension>& coords,
      std::array<int_t,dimension>& bits)
  {
    if(bits[1] == 0){
      if(bits[0] == 1){
        coords[0] = n-1 - coords[0];
        coords[1] = n-1 - coords[1];
      }
      // Swap X-Y or Z
      int t = coords[0];
      coords[0] = coords[1];
      coords[1] = t;
    }
  }


  /* 2D Hilbert key
   * Hilbert key is always generated to the max_depth and then truncated
   * */
  template<
    typename S>
  hilbert_id(
    const std::array<point__<S, dimension>, 2>& range,
    const point__<S, dimension>& p,
    size_t depth)
  : id_(int_t(1) << max_depth * dimension + bits%dimension)
  {
    assert(depth <= max_depth);

    std::array<int_t, dimension> coords{};
    //std::fill(coords.begin(),coords.end(),0);

    // Convert the position to integer
    for(size_t i = 0; i < dimension; ++i)
    {
      S min = range[0][i];
      S scale = range[1][i] - min;
      coords[i] = (p[i] - min)/scale * (int_t(1) << ((bits - 1)/dimension));
    }

    if(dimension == 1)
    {
      std::cout<<coords[0]<<std::endl;
      assert(id_ & 1UL<<max_depth);
      id_ |= coords[0]>>dimension;
      id_ >>= (max_depth-depth);
      std::cout<<"k: "<<std::bitset<64>(id_)<<std::endl<<std::flush;
      return;
    }

    //std::cout<<"b: "<<std::bitset<64>(id_)<<std::endl;
    //std::cout<<"max_depth="<<max_depth<<" depth="<<depth<<std::endl;
    int_t mask = static_cast<int_t>(1)<<(max_depth);
    for(int_t s = mask>>1 ; s > 0; s >>= 1)
    {
      std::array<int_t,dimension> bits;
      for(size_t j = 0 ; j < dimension; ++j)
      {
        bits[j] = (s & coords[j]) > 0;
      }
      id_ += s * s * ((3*bits[0]) ^ bits[1]);
      rotation(s,coords,bits);
    }
    // Then truncate the key to the depth
    id_ >>= (max_depth-depth)*dimension;
    //std::cout<<"pos: "<<p<<std::endl;
    //std::cout<<"k: "<<std::bitset<64>(id_)<<std::endl;

    // Be sure the head bit is one
    //assert((id_ & (int_t(1)<<depth*dimension+bits%dimension)) > 0);
  }

  /*!
    Construct the id from an array of dimensions and range for each
    dimension. Construct with the max depth available.
   */
  template<
    typename S
  >
  hilbert_id(
    const std::array<point__<S, dimension>, 2>& range,
    const point__<S, dimension>& p)
  : hilbert_id(range,p,max_depth)
  {}


  constexpr hilbert_id(const hilbert_id& bid)
  : id_(bid.id_)
  {}

  /*!
    Get the root id (depth 0).
   */
  static
  constexpr
  hilbert_id
  root()
  {
    return hilbert_id(int_t(1));
  }

  /*!
    Get the null id.
   */
  static
  constexpr
  hilbert_id
  null()
  {
    return hilbert_id(0);
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

  hilbert_id&
  operator=(
    const hilbert_id& bid
  )
  {
    id_ = bid.id_;
    return *this;
  }

  constexpr
  bool
  operator==(
    const hilbert_id& bid
  ) const
  {
    return id_ == bid.id_;
  }

  constexpr
  bool
  operator!=(
    const hilbert_id& bid
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
  hilbert_id
  parent() const
  {
    return hilbert_id(id_ >> dimension);
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
    ostr<<std::bitset<64>(id_);
  }

  int_t
  value_() const
  {
    return id_;
  }

  bool
  operator<(
    const hilbert_id& bid
  ) const
  {
    return id_ < bid.id_;
  }

  /*!
    Convert this id to coordinates in range.
   */
/*  template<
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
  }*/

  /**
   * @brief Compute the range of a branch from its key
   * The space is recursively decomposed regarding the dimension
   */
/*  template<
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
    int_t root = hilbert_id::root().id_;

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
  }*/

private:
  int_t id_;

  constexpr
  hilbert_id(
    int_t id
  )
  : id_(id)
  {}
};

// output for hilbert id
template<
  typename T,
  size_t D>
std::ostream&
operator<<(
  std::ostream& ostr,
  const hilbert_id<T,D>& k)
{
  k.output_(ostr);
  return ostr;
}

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_hilbert_id_h