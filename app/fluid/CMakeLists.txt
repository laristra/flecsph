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

#------------------------------------------------------------------------------#
# Includes
#------------------------------------------------------------------------------#
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/app/fluid/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(
  ${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#
# Fluid
add_executable(fluid 
  main.cc 
  main_driver.cc 
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

# Fluid generator
add_executable(fluid_generator 
  generator/main.cc 
)

target_link_libraries(fluid_generator ${MPI_LIBRARIES} 
  ${HDF5_LIBRARIES} ${H5hut_LIBRARIES})
target_link_libraries(fluid ${FleCSPH_LIBRARIES})

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)
