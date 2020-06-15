/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

/*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file tree.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Description of the tree policy used in the tree_topology from FleCSI.
 */

#ifndef tree_h
#define tree_h

#include <vector>

//#warning "CHANGE TO FLECSI ONE"
#include "space_vector.h"
#include "tree_topology/filling_curve.h"
#include "tree_topology/tree_topology.h"
//#include "utils.h"

#include "body.h"
#include "node.h"
#include <boost/multiprecision/cpp_int.hpp>

using namespace flecsi;
using boost::multiprecision::uint128_t;

#ifdef KEY_INTEGER_TYPE
using key_type_t = KEY_INTEGER_TYPE;
#else
using key_type_t = uint64_t;
// using key_type_t = uint128_t;
#endif

namespace flecsi {
namespace execution {
void specialization_driver(int argc, char * argv[]);
void driver(int argc, char * argv[]);
} // namespace execution
} // namespace flecsi

class tree_policy
{
public:
  using tree_t = flecsi::topology::tree_topology<tree_policy>;
  using key_int_t = key_type_t;
  static const size_t dimension = gdimension;
  using element_t = type_t;
  using key_t = flecsi::morton_curve_u<dimension, key_type_t>;
  using point_t = flecsi::space_vector_u<element_t, dimension>;
  using geometry_t = flecsi::topology::tree_geometry<element_t, gdimension>;
  using entity_t = body_u<key_t>;
  using cofm_t = node_u<key_t,1>; 
}; // class tree_policy

using tree_topology_t = flecsi::topology::tree_topology<tree_policy>;
using tree_geometry_t = flecsi::topology::tree_geometry<type_t, gdimension>;
using point_t = tree_topology_t::point_t;
using node = tree_topology_t::cofm_t;
using key_type = tree_topology_t::key_t;
using body = tree_topology_t::entity_t;

using range_t = std::array<point_t, 2>;

/* TODO: do we still need these?
inline
bool
operator<(
    const point_t& p,
    const point_t& q)
{
  for(size_t i=0;i<gdimension;++i)
    if(p[i]>q[i])
      return false;
  return true;
}

inline
bool
operator>(
    const point_t& p,
    const point_t& q)
{
  for(size_t i=0;i<gdimension;++i)
    if(p[i]<q[i])
      return false;
  return true;
}

inline
point_t
operator*(
    const point_t& p,
    const point_t& q)
{
  point_t r = p;
  for(size_t i=0;i<gdimension;++i)
    r[i] *= q[i];
  return r;
}
*/

#endif // tree_h
