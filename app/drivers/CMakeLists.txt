#----------------------------------------------------------------------------#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#
option(ENABLE_MPI_TESTS "Enable unit testing with MPI" ON)

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
include_directories(${CMAKE_SOURCE_DIR}/include/physics/eos)
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
install(TARGETS hydro_1d RUNTIME DESTINATION bin)

# MPILEGION version: 2D
add_executable(hydro_2d
  hydro/main.cc
  hydro/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)


target_link_libraries(hydro_2d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_2d PUBLIC -DEXT_GDIMENSION=2)
install(TARGETS hydro_2d RUNTIME DESTINATION bin)

# MPILEGION version: 3D
add_executable(hydro_3d
  hydro/main.cc
  hydro/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(hydro_3d ${FleCSPH_LIBRARIES})
target_compile_definitions(hydro_3d PUBLIC -DEXT_GDIMENSION=3)
install(TARGETS hydro_3d RUNTIME DESTINATION bin)


#------------------------------------------------------------------------------#
# WVT drivers
#------------------------------------------------------------------------------#

# 2D-version
add_executable(wvt_2d
  wvt/main.cc
  wvt/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)
target_link_libraries(wvt_2d ${FleCSPH_LIBRARIES})
target_compile_definitions(wvt_2d PUBLIC -DEXT_GDIMENSION=2)
install(TARGETS wvt_2d RUNTIME DESTINATION bin/drivers)


# 3D-version
add_executable(wvt_3d
  wvt/main.cc
  wvt/main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)
target_link_libraries(wvt_3d ${FleCSPH_LIBRARIES})
target_compile_definitions(wvt_3d PUBLIC -DEXT_GDIMENSION=3)
install(TARGETS wvt_3d RUNTIME DESTINATION bin/drivers)


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
install(TARGETS newtonian_3d RUNTIME DESTINATION bin)

#------------------------------------------------------------------------------#
# collapse test, call the default parameter file
#------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
cinch_add_unit(collapse_test
  SOURCES
    test/collapse.cc
    newtonian/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/collapse_nx10.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=3
  POLICY MPI
)

add_custom_target(
    collapse_3d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/collapse_3d_generator ${PROJECT_SOURCE_DIR}/data/collapse_nx10.par
    COMMAND mv collapse_nx10.h5part ../../test/drivers
    DEPENDS ../id_generators/collapse_3d_generator
)
add_dependencies(collapse_test collapse_3d_generator_test)

if(ENABLE_MPI_TESTS)
cinch_add_unit(collapse_MPI_test
  SOURCES
    test/collapse.cc
    newtonian/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/collapse_nx10.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=3
  POLICY MPI
   THREADS 4
)

add_dependencies(collapse_MPI_test collapse_3d_generator_test)
endif()  # ENABLE_MPI_TESTS
endif()  # ENABLE_UNIT_TESTS


#------------------------------------------------------------------------------#
# sodtube test, call the default parameter file
#------------------------------------------------------------------------------#
if(ENABLE_UNIT_TESTS)

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
)

add_custom_target(
    sodtube_1d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/sodtube_1d_generator ${PROJECT_SOURCE_DIR}/data/sodtube_t1_n100.par
    COMMAND mv sodtube_t1_n100.h5part ../../test/drivers
    DEPENDS ../id_generators/sodtube_1d_generator
)
add_dependencies(sodtube_test sodtube_1d_generator_test)

if(ENABLE_MPI_TESTS)
cinch_add_unit(sodtube_MPI_test
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
   THREADS 4
)

add_dependencies(sodtube_MPI_test sodtube_1d_generator_test)
endif()  # ENABLE_MPI_TESTS
endif()  # ENABLE_UNIT_TESTS
#------------------------------------------------------------------------------#
# RT test, call the default parameter file
#------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)

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

add_custom_target(
    RT_2d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/RT_2d_generator ${PROJECT_SOURCE_DIR}/data/RT_2d.par
    COMMAND mv RT_n30.h5part ../../test/drivers
    DEPENDS ../id_generators/RT_2d_generator
)
add_dependencies(RT_test RT_2d_generator_test)
endif()

#------------------------------------------------------------------------------#
# KH test: uses KH_2d.par from data directory
#------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)

cinch_add_unit(KH_test
  SOURCES
    test/KH.cc
    hydro/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/KH_2d.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=2
  POLICY MPI
  #THREADS 2
)

add_custom_target(
    KH_2d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/KH_2d_generator ${PROJECT_SOURCE_DIR}/data/KH_2d.par
    COMMAND mv KH_t1_n30.h5part ../../test/drivers
    DEPENDS ../id_generators/KH_2d_generator
)
add_dependencies(KH_test KH_2d_generator_test)

endif() # ENABLE_UNIT_TESTS

#------------------------------------------------------------------------------#
# sedov test with the default parameter file
#------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
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

add_custom_target(
    sedov_2d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/sedov_2d_generator ${PROJECT_SOURCE_DIR}/data/sedov_nx20.par
    COMMAND mv sedov_nx20.h5part ../../test/drivers
    DEPENDS ../id_generators/sedov_2d_generator
)
add_dependencies(sedov_test sedov_2d_generator_test)
endif()

#------------------------------------------------------------------------------#
# relaxation test for the "mesa" potential in 3D using mesa_nx20.par file
#------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
cinch_add_unit(mesa_relaxation_test
  SOURCES
    test/mesa.cc
    hydro/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/mesa_nx20.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
  DEFINES
    -DEXT_GDIMENSION=3
  POLICY MPI
  #THREADS 2
)

add_custom_target(
    mesa_3d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/sedov_3d_generator ${PROJECT_SOURCE_DIR}/data/mesa_nx20.par
    COMMAND mv mesa_nx20.h5part ../../test/drivers
    DEPENDS ../id_generators/sedov_3d_generator
)
add_dependencies(mesa_relaxation_test mesa_3d_generator_test)
endif()

#------------------------------------------------------------------------------#
# noh test, with the default parameter file
#------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
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

add_custom_target(
    noh_2d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/noh_2d_generator ${PROJECT_SOURCE_DIR}/data/noh_nx20.par
    COMMAND mv noh_nx20.h5part ../../test/drivers
    DEPENDS ../id_generators/noh_2d_generator
)
add_dependencies(noh_test noh_2d_generator_test)
endif()

#------------------------------------------------------------------------------#
# implosion test
#------------------------------------------------------------------------------#

if(ENABLE_UNIT_TESTS)
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

add_custom_target(
    implosion_2d_generator_test
    COMMAND ${MPIEXEC} -n 1 ../id_generators/implosion_2d_generator ${PROJECT_SOURCE_DIR}/data/implosion_nx20.par
    COMMAND mv implosion_nx20.h5part ../../test/drivers
    DEPENDS ../id_generators/implosion_2d_generator
)
add_dependencies(implosion_test implosion_2d_generator_test)
endif()
