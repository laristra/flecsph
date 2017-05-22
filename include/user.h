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
 * @file user.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Main description to switch from SodTube to BNS  
 */

#ifndef user_h
#define user_h

#define BNS
//#define SODTUBE
#define OUTPUT
//#define OUTPUTGRAPH

#ifdef BNS
static const size_t gdimension = 3;
using type_t = double;
#endif

#ifdef SODTUBE
static const size_t gdimension = 1;
using type_t = double; 
#endif

#endif // user_h

