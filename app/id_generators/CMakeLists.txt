#----------------------------------------------------------------------------#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Debug and release flags
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG
  "-DPARALLEL_IO -Wno-sign-compare -Wno-reorder -Wno-narrowing -Wno-deprecated-declaration -ftree-vectorize -ffast-math -msse2 -g -O2 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")
set(CMAKE_CXX_FLAGS_RELEASE
  "-DPARALLEL_IO -Wno-sign-compare -Wno-reorder -Wno-narrowing -Wno-deprecated-declaration -ftree-vectorize -ffast-math -msse2 -O3 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")

# includes

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/physics)
include_directories(${CMAKE_SOURCE_DIR}/app/drivers/include)
include_directories(${CMAKE_SOURCE_DIR}/app/id_generators/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)

#------------------------------------------------------------------------------#
# Executables
#------------------------------------------------------------------------------#

add_executable(sodtube_1d_generator     sodtube/main.cc)
add_executable(sodtube_2d_generator     sodtube/main.cc)
add_executable(sodtube_3d_generator     sodtube/main.cc)
add_executable(sedov_2d_generator       sedov/main.cc)
add_executable(sedov_3d_generator       sedov/main.cc)
add_executable(noh_2d_generator         noh/main.cc)
add_executable(noh_3d_generator         noh/main.cc)
add_executable(implosion_2d_generator   implosion/main.cc)
add_executable(implosion_3d_generator   implosion/main.cc)
add_executable(collapse_3d_generator    collapse/main.cc)
add_executable(wtunnel_2d_generator     wtunnel/main.cc)
add_executable(wtunnel_3d_generator     wtunnel/main.cc)
add_executable(KH_2d_generator          KH/main.cc)
add_executable(KH_3d_generator          KH/main.cc)
add_executable(RT_2d_generator          RT/main.cc)
add_executable(RT_3d_generator          RT/main.cc)
add_executable(fluid_2d_generator          fluid/main.cc)
add_executable(fluid_3d_generator          fluid/main.cc)


# Install
install(TARGETS
    sodtube_1d_generator
    sodtube_2d_generator
    sodtube_3d_generator
    sedov_2d_generator
    sedov_3d_generator
    noh_2d_generator
    noh_3d_generator
    implosion_2d_generator
    implosion_3d_generator
    wtunnel_2d_generator
    wtunnel_3d_generator
    KH_2d_generator
    KH_3d_generator
    RT_2d_generator
    RT_3d_generator
    fluid_2d_generator
    fluid_3d_generator
    RUNTIME
    DESTINATION bin/id_generators)

target_link_libraries(sodtube_1d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sodtube_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sodtube_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sedov_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sedov_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(noh_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(noh_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(implosion_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(implosion_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(collapse_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(wtunnel_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(wtunnel_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(KH_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(KH_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(RT_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(RT_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(fluid_2d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(fluid_3d_generator
  ${MPI_LIBRARIES} ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})

target_compile_definitions(sodtube_1d_generator     PUBLIC -DEXT_GDIMENSION=1)
target_compile_definitions(sodtube_2d_generator     PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(sodtube_3d_generator     PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(sedov_2d_generator       PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(sedov_3d_generator       PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(noh_2d_generator         PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(noh_3d_generator         PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(implosion_2d_generator   PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(implosion_3d_generator   PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(collapse_3d_generator    PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(wtunnel_2d_generator     PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(wtunnel_3d_generator     PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(KH_2d_generator          PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(KH_3d_generator          PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(RT_2d_generator          PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(RT_3d_generator          PUBLIC -DEXT_GDIMENSION=3)
target_compile_definitions(fluid_2d_generator          PUBLIC -DEXT_GDIMENSION=2)
target_compile_definitions(fluid_3d_generator          PUBLIC -DEXT_GDIMENSION=3)
