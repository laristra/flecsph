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

include_directories(${CMAKE_SOURCE_DIR}/app/sedov/include)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/bodies)
include_directories(${CMAKE_SOURCE_DIR}/include/physics)
include_directories(${CMAKE_SOURCE_DIR}/include/physics/default)
include_directories(${CMAKE_SOURCE_DIR}/mpisph/test)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)


#------------------------------------------------------------------------------#
# tree topology tests
#------------------------------------------------------------------------------#

cinch_add_unit(sedov_test
  SOURCES
    sedov.cc
    ${CMAKE_SOURCE_DIR}/app/sedov/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc  
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/hdf5_sedov.h5part
  LIBRARIES 
    ${FleCSPH_LIBRARIES}
)

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
