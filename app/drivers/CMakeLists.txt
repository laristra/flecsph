#----------------------------------------------------------------------------#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Debug and release flags
#------------------------------------------------------------------------------#

# For vectorization: -fopt-info-vec-missed/optimized/note/all
set(CMAKE_CXX_FLAGS_DEBUG
  "-mavx2 -msse4 -ftree-vectorize -march=native \
  -pthread -lpthread -DPARALLEL_IO -Wno-sign-compare -Wno-reorder -Wno-narrowing \
  -Wno-deprecated-declarations  -ffast-math -g -O3 -Wall -Wno-return-type \
-Wno-unused -Wno-comment -Wno-parentheses -lm")
set(CMAKE_CXX_FLAGS_RELEASE
 "-march=native -pthread -lpthread -DPARALLEL_IO -Wno-sign-compare -Wno-reorder \
  -Wno-narrowing -Wno-deprecated-declarations -ftree-vectorize -ffast-math -msse2 \
  -O3 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")

# includes

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/physics)
include_directories(${CMAKE_SOURCE_DIR}/include/physics/specific)
include_directories(${CMAKE_SOURCE_DIR}/app/drivers/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)

#------------------------------------------------------------------------------#
# Hydro drivers without gravity
#------------------------------------------------------------------------------#

# MPILEGION version: 1D
add_executable(hydro_1d
  hydro/main.cc
  hydro/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)


target_link_libraries(hydro_1d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_1d PUBLIC -DEXT_GDIMENSION=1)
install(TARGETS hydro_1d RUNTIME DESTINATION bin/drivers)

# MPILEGION version: 2D
add_executable(hydro_2d
  hydro/main.cc
  hydro/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)


target_link_libraries(hydro_2d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_2d PUBLIC -DEXT_GDIMENSION=2)
install(TARGETS hydro_2d RUNTIME DESTINATION bin/drivers)

# MPILEGION version: 3D
add_executable(hydro_3d
  hydro/main.cc
  hydro/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(hydro_3d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_3d PUBLIC -DEXT_GDIMENSION=3)
install(TARGETS hydro_3d RUNTIME DESTINATION bin/drivers)

#------------------------------------------------------------------------------#
# Hydro drivers with Newtonian gravity
#------------------------------------------------------------------------------#

# 3D-version
add_executable(newtonian_3d
  newtonian/main.cc
  newtonian/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)
target_link_libraries(newtonian_3d ${FleCSPH_LIBRARIES})
#### THESE ARE NOT BUILT PROPERLY as of now -- commenting them out (2018-10-05, O.K.)
# target_link_libraries(newtonian_3d ${CMAKE_SOURCE_DIR}/third-party-libraries/stellar_collapse/libTABEOS_RF.a)
# target_link_libraries(newtonian_3d ${CMAKE_SOURCE_DIR}/third-party-libraries/stellar_collapse/libTABEOS_SC.a)
target_compile_definitions(newtonian_3d PUBLIC -DEXT_GDIMENSION=3)
install(TARGETS newtonian_3d RUNTIME DESTINATION bin/drivers)


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
  POLICY MPI
  # THREADS 2
)

if(ENABLE_UNIT_TESTS)

add_custom_target(
    sodtube_1d_generator_test
    COMMAND mpirun -n 1 ../id_generators/sodtube_1d_generator ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par
    COMMAND mv sodtube_t1_n100.h5part ../../test/drivers
    DEPENDS ../id_generators/sodtube_1d_generator
)
add_dependencies(sodtube_test sodtube_1d_generator_test)
endif()

#------------------------------------------------------------------------------#
# RT test, call the default parameter file
#------------------------------------------------------------------------------#

cinch_add_unit(RT_test
  SOURCES
    test/RT.cc
    hydro/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/RT_2d.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=2
  POLICY MPI
  #THREADS 2
)

if(ENABLE_UNIT_TESTS)

add_custom_target(
    RT_2d_generator_test
    COMMAND mpirun -n 1 ../id_generators/RT_2d_generator ${PROJECT_SOURCE_DIR}/data/RT_2d.par
    COMMAND mv RT_n30.h5part ../../test/drivers
    DEPENDS ../id_generators/RT_2d_generator
)
add_dependencies(RT_test RT_2d_generator_test)
endif()

#------------------------------------------------------------------------------#
# sedov test with the default parameter file
#------------------------------------------------------------------------------#

cinch_add_unit(sedov_test
  SOURCES
    test/sedov.cc
    hydro/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/sedov_nx20.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=2
  POLICY MPI
  #THREADS 2
)

if(ENABLE_UNIT_TESTS)
add_custom_target(
    sedov_2d_generator_test
    COMMAND mpirun -n 1 ../id_generators/sedov_2d_generator ${PROJECT_SOURCE_DIR}/data/sedov_nx20.par
    COMMAND mv sedov_nx20.h5part ../../test/drivers
    DEPENDS ../id_generators/sedov_2d_generator
)
add_dependencies(sedov_test sedov_2d_generator_test)
endif()

#------------------------------------------------------------------------------#
# noh test, with the default parameter file
#------------------------------------------------------------------------------#

cinch_add_unit(noh_test
  SOURCES
    test/noh.cc
    hydro/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/noh_nx20.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=2
  POLICY MPI
  #THREADS 2
)

if(ENABLE_UNIT_TESTS)
add_custom_target(
    noh_2d_generator_test
    COMMAND mpirun -n 1 ../id_generators/noh_2d_generator ${PROJECT_SOURCE_DIR}/data/noh_nx20.par
    COMMAND mv noh_nx20.h5part ../../test/drivers
    DEPENDS ../id_generators/noh_2d_generator
)
add_dependencies(noh_test noh_2d_generator_test)
endif()

#------------------------------------------------------------------------------#
# implosion test
#------------------------------------------------------------------------------#

cinch_add_unit(implosion_test
  SOURCES
    test/implosion.cc
    hydro/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/implosion_nx20.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=2
  POLICY MPI
  #THREADS 2
)

if(ENABLE_UNIT_TESTS)
add_custom_target(
    implosion_2d_generator_test
    COMMAND mpirun -n 1 ../id_generators/implosion_2d_generator ${PROJECT_SOURCE_DIR}/data/implosion_nx20.par
    COMMAND mv implosion_nx20.h5part ../../test/drivers
    DEPENDS ../id_generators/implosion_2d_generator
)
add_dependencies(implosion_test implosion_2d_generator_test)
endif()
