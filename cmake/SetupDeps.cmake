#------------------------------------------------------------------------------#
# Find GoogleTest
#------------------------------------------------------------------------------#
if(ENABLE_UNIT_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
endif()
#------------------------------------------------------------------------------#
# Find Threads
#------------------------------------------------------------------------------#
find_package(Threads REQUIRED)

#------------------------------------------------------------------------------#
# Find OpenMP
#------------------------------------------------------------------------------#
find_package(OpenMP REQUIRED COMPONENTS CXX)
# list (APPEND FleCSPH_LIBRARIES OpenMP::OpenMP_CXX)

#------------------------------------------------------------------------------#
# MPI
#------------------------------------------------------------------------------#
if(DEFINED ENV{CRAY_MPICH2_DIR})
  find_library(MPICH_LIB
               NAMES mpich
               HINTS $ENV{CRAY_MPICH2_DIR}/lib)
  set(MPI_mpich_LIBRARY ${MPICH_LIB})
  set(MPI_CXX_LIB_NAMES "mpich" CACHE STRING "" FORCE)
  set(MPI_CXX_HEADER_DIR $ENV{CRAY_MPICH2_DIR}/include)
endif()
find_package(MPI REQUIRED COMPONENTS CXX)
# list (APPEND FleCSPH_LIBRARIES MPI::MPI_CXX )

#------------------------------------------------------------------------------#
# Add GSL
#------------------------------------------------------------------------------#
find_package(GSL REQUIRED)

#------------------------------------------------------------------------------#
# Add Boost
#------------------------------------------------------------------------------#
find_package(Boost REQUIRED)
# include_directories(${Boost_INCLUDE_DIR})

#------------------------------------------------------------------------------#
# HDF5
#------------------------------------------------------------------------------#
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)

#------------------------------------------------------------------------------#
# FleCSI
#------------------------------------------------------------------------------#
find_package(FleCSI REQUIRED)
