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
include_directories(${CMAKE_SOURCE_DIR}/app/drivers/include)
include_directories(${CMAKE_SOURCE_DIR}/app/id_generators/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# MPILEGION version: 1D
add_executable(hydro_1d
  hydro/main.cc 
  hydro/main_driver.cc 
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(hydro_1d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_1d PUBLIC -DEXT_GDIMENSION=1)

# MPILEGION version: 2D
add_executable(hydro_2d
  hydro/main.cc 
  hydro/main_driver.cc 
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(hydro_2d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_2d PUBLIC -DEXT_GDIMENSION=2)

# MPILEGION version: 3D
add_executable(hydro_3d
  hydro/main.cc 
  hydro/main_driver.cc 
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(hydro_3d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_3d PUBLIC -DEXT_GDIMENSION=3)

#------------------------------------------------------------------------------#
# sodtube test, call the default parameter file 
#------------------------------------------------------------------------------#

cinch_add_unit(sodtube_test
  SOURCES
    test/sodtube.cc
    hydro/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc 
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par 
  LIBRARIES 
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=1
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
