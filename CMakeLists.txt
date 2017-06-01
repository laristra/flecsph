#------------------------------------------------------------------------------#
# Require some version of cmake
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 2.8)
project(tree CXX)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

if(NOT "$ENV{CINCH_SOURCE_DIR}" STREQUAL "")
  set(CINCH_SOURCE_DIR $ENV{CINCH_SOURCE_DIR})
elseif(EXISTS ${CMAKE_SOURCE_DIR}/cinch)
  set(CINCH_SOURCE_DIR ${CMAKE_SOURCE_DIR}/cinch)
endif()



set(FLECSI_RUNTIME_MODEL "serial" CACHE STRING "The flecsi runtime")

include(${CINCH_SOURCE_DIR}/cmake/ProjectLists.txt)
