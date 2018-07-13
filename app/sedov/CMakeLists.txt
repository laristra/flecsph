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
include_directories(${CMAKE_SOURCE_DIR}/app/sedov/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)
#include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/local/include)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# MPILEGION version
add_executable(sedov
  main.cc 
  main_driver.cc 
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

# MPILEGION version
add_executable(sedov_generator 
  generator/main.cc 
)

target_link_libraries(sedov_generator ${MPI_LIBRARIES} 
  ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sedov ${FleCSPH_LIBRARIES})

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)

#------------------------------------------------------------------------------#
# sedov test, call the default parameter file 
#------------------------------------------------------------------------------#

cinch_add_unit(sedov_test
  SOURCES
    test/sedov.cc
    ${CMAKE_SOURCE_DIR}/app/sedov/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc 
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/sedov_sqn100.par 
  LIBRARIES 
    ${FleCSPH_LIBRARIES}
)

if(ENABLE_UNIT_TESTS)
# run the generator with default file 
add_custom_target(
    sedov_generator_test 
    COMMAND sedov_generator ${PROJECT_SOURCE_DIR}/data/sedov_sqn100.par
    COMMAND mv sedov_sqn100.h5part ../../test/sedov
)
add_dependencies(sedov_test sedov_generator_test)
endif()