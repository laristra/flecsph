#~----------------------------------------------------------------------------~#
#  Copyright (c) 2017 Los Alamos National Security, LLC
#  All rights reserved.
#~----------------------------------------------------------------------------~#

#~----------------------------------------------------------------------------~#
# 
# /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
# /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
# /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
# /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
# /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
# /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
# /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
# //       ///  //////   //////  ////////  //       //      //  
#
#~---------------------------------------------------------------------------~*#

include_directories(${CMAKE_SOURCE_DIR}/app/sodtube/include)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/bodies)
include_directories(${CMAKE_SOURCE_DIR}/include/physics)
include_directories(${CMAKE_SOURCE_DIR}/include/physics/default)
include_directories(${CMAKE_SOURCE_DIR}/mpisph/test)
include_directories(${CMAKE_SOURCE_DIR}/mpisph/)



#------------------------------------------------------------------------------#
# tree topology tests
#------------------------------------------------------------------------------#

cinch_add_unit(sodtube_test
  SOURCES
    sodtube.cc
    ${CMAKE_SOURCE_DIR}/app/sodtube/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc 
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.h5part 
    ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par 
  LIBRARIES 
    ${FleCSPH_LIBRARIES}
)

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
