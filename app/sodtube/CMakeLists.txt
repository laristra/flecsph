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
include_directories(${CMAKE_SOURCE_DIR}/app/sodtube/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# MPILEGION version
add_executable(sodtube 
  main.cc 
  main_driver.cc 
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

# MPILEGION version
add_executable(sodtube_generator 
  generator/main.cc 
)

target_link_libraries(sodtube_generator ${MPI_LIBRARIES} 
  ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sodtube ${FleCSPH_LIBRARIES})

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)

#------------------------------------------------------------------------------#
# sodtube test, call the default parameter file 
#------------------------------------------------------------------------------#

cinch_add_unit(sodtube_test
  SOURCES
    test/sodtube.cc
    ${CMAKE_SOURCE_DIR}/app/sodtube/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc 
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par 
  LIBRARIES 
    ${FleCSPH_LIBRARIES}
)

if(ENABLE_UNIT_TESTS)
# run the generator with default file 
add_custom_target(
    sodtube_generator_test 
    COMMAND sodtube_generator ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par
    COMMAND mv sodtube_t1_n100.h5part ../../test/sodtube
)
add_dependencies(sodtube_test sodtube_generator_test)
endif()