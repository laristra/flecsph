#----------------------------------------------------------------------------#
# Copyright (c) 2017 Los Alamos National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Debug and release flags     
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG 
  "-g -O2 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")
set(CMAKE_CXX_FLAGS_RELEASE 
  "-O3 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses") 

# includes 

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/physics)
include_directories(${CMAKE_SOURCE_DIR}/app/id_generators/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

add_executable(sodtube_generator sodtube/main.cc)
add_executable(sedov_generator   sedov/main.cc)
add_executable(noh_generator     noh/main.cc)

target_link_libraries(sodtube_generator 
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sedov_generator 
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(noh_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})

