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
# Find FleCSI
#------------------------------------------------------------------------------#

find_package(FleCSI REQUIRED)
include_directories(${FleCSI_INCLUDE_DIR})
set(FleCSPH_LIBRARIES ${FleCSI_LIBRARIES})

#------------------------------------------------------------------------------#
# Pthreads for parallelism in the tree
#------------------------------------------------------------------------------#

#find_package(Threads REQUIRED)
#include_directories(${THREADS_INCLUDE_PATH})

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
# HDF5
#------------------------------------------------------------------------------#

find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIR})
list(APPEND FleCSPH_LIBRARIES ${HDF5_LIBRARIES})

#------------------------------------------------------------------------------#
# Add OpenMP
#------------------------------------------------------------------------------#

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
