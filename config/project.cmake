cinch_minimum_required(1.0)

project(tree)

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
# Find FleCSI
#------------------------------------------------------------------------------#

find_package(FleCSI REQUIRED)
include_directories(${FleCSI_INCLUDE_DIR})
message(STATUS ${FleCSI_INCLUDE_DIR})
set(FleCSPH_LIBRARIES ${FleCSI_LIBRARIES})

#------------------------------------------------------------------------------#
# Pthreads for parallelism in the tree
#------------------------------------------------------------------------------#

#find_package(Threads REQUIRED)
#include_directories(${THREADS_INCLUDE_PATH})

#------------------------------------------------------------------------------#
# HDF5
#------------------------------------------------------------------------------#

#find_package(HDF5)
set(HDF5_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party-libraries/local/include")
set(HDF5_LIBRARIES "${CMAKE_SOURCE_DIR}/third-party-libraries/local/lib/libhdf5.so")

include_directories(${HDF5_INCLUDE_DIR})
list(APPEND FleCSPH_LIBRARIES ${HDF5_LIBRARIES})

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
set(H5HUT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party-libraries/local/include")
set(H5HUT_LIBRARIES "${CMAKE_SOURCE_DIR}/third-party-libraries/local/lib/libH5hut.so")

include_directories(${H5HUT_INCLUDE_DIR})
list(APPEND FleCSPH_LIBRARIES ${H5HUT_LIBRARIES})

#------------------------------------------------------------------------------#
# Add library target
#------------------------------------------------------------------------------#

cinch_add_library_target(mpisph "mpisph")

#------------------------------------------------------------------------------#
# Add application targets
#------------------------------------------------------------------------------#

cinch_add_application_directory("app/fluid")
cinch_add_application_directory("app/sodtube")
