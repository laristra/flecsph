#------------------------------------------------------------------------------#
# The cinch version 
#------------------------------------------------------------------------------#
cinch_minimum_required(1.0)

#------------------------------------------------------------------------------#
# The project name 
#------------------------------------------------------------------------------#
project(flecsph)


#------------------------------------------------------------------------------#
# Set the minimum Cinch version
#------------------------------------------------------------------------------#
set(CINCH_HEADER_SUFFIXES "\\.h")

#------------------------------------------------------------------------------#
# If a C++14 compiler is available, then set the appropriate flags
#------------------------------------------------------------------------------#
include(cxx14)
check_for_cxx14_compiler(CXX14_COMPILER)
if(CXX14_COMPILER)
  enable_cxx14()
else()
  message(FATAL_ERROR "C++14 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# Loading cinch extras
#------------------------------------------------------------------------------#
cinch_load_extras()


#------------------------------------------------------------------------------#
# Find FleCSI
#------------------------------------------------------------------------------#
find_package(FleCSI REQUIRED)
include_directories(${FleCSI_INCLUDE_DIR})
message(STATUS ${FleCSI_INCLUDE_DIR})
set(FleCSPH_LIBRARIES ${FleCSI_LIBRARIES})
message(STATUS ${FleCSPH_LIBRARIES})

#------------------------------------------------------------------------------#
# MPI
#------------------------------------------------------------------------------#
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
list(APPEND FleCSPH_LIBRARIES ${MPI_LIBRARIES})

#------------------------------------------------------------------------------#
# Legion
#------------------------------------------------------------------------------#
find_package(Legion REQUIRED)
include_directories(${Legion_INCLUDE_DIRS})
list(APPEND FleCSPH_LIBRARIES ${Legion_LIBRARY} ${REALM_LIBRARY})

#------------------------------------------------------------------------------#
# Pthreads for parallelism in the tree
#------------------------------------------------------------------------------#
#find_package(Threads REQUIRED)
#include_directories(${THREADS_INCLUDE_PATH})

#------------------------------------------------------------------------------#
# Add OpenMP
#------------------------------------------------------------------------------#
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#------------------------------------------------------------------------------#
# H5HUT
#------------------------------------------------------------------------------#
# Set by hand right now
find_package(H5hut REQUIRED)
#set(HDF5_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party-libraries/local/include")
#set(HDF5_LIBRARIES "${CMAKE_SOURCE_DIR}/third-party-libraries/local/lib/libH5hut.so")
message(STATUS ${H5hut_LIBRARIES})
message(STATUS ${H5hut_INCLUDE_DIRS})
list(APPEND FleCSPH_LIBRARIES ${H5hut_LIBRARIES})
include_directories(${H5hut_INCLUDE_DIRS})

#------------------------------------------------------------------------------#
# HDF5
#------------------------------------------------------------------------------#
#for static libH4hut
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
#set(HDF5_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party-libraries/local/include")
#set(HDF5_LIBRARIES "${CMAKE_SOURCE_DIR}/third-party-libraries/local/lib/libhdf5.so")
message(STATUS ${HDF5_LIBRARIES})
message(STATUS ${HDF5_INCLUDE_DIR})
list(APPEND FleCSPH_LIBRARIES ${HDF5_LIBRARIES})
include_directories(${HDF5_INCLUDE_DIR})

#------------------------------------------------------------------------------#
# HDF5 ScalingFramework
#------------------------------------------------------------------------------#
set(HSF_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests")
message(STATUS ${HSF_INCLUDE_DIR})
include_directories(${HSF_INCLUDE_DIR})



#------------------------------------------------------------------------------#
# Add mpisph tests
#------------------------------------------------------------------------------#
cinch_add_application_directory("mpisph/")

#------------------------------------------------------------------------------#
# Add application targets
#------------------------------------------------------------------------------#
cinch_add_application_directory("app/bns_3D")
cinch_add_application_directory("app/sedov")
cinch_add_application_directory("app/noh")
cinch_add_application_directory("app/fluid_3D")
cinch_add_application_directory("app/fluid_2D")
cinch_add_application_directory("app/sodtube")
cinch_add_application_directory("app/bwd")
