#------------------------------------------------------------------------------#
# The cinch version
#------------------------------------------------------------------------------#
cinch_minimum_required(VERSION v1.0)

#------------------------------------------------------------------------------#
# The project name
#------------------------------------------------------------------------------#
project(flecsph)


#------------------------------------------------------------------------------#
# Set the minimum Cinch version
#------------------------------------------------------------------------------#
set(CINCH_HEADER_SUFFIXES "\\.h")

#------------------------------------------------------------------------------#
# A C++17-capable compiler is required
#------------------------------------------------------------------------------#
include(cxx17)
check_for_cxx17_compiler(CXX17_COMPILER)
if(CXX17_COMPILER)
  enable_cxx17()
else()
  message(FATAL_ERROR "C++17 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# Loading cinch extras
#------------------------------------------------------------------------------#
cinch_load_extras()

#------------------------------------------------------------------------------#
# Legion
#------------------------------------------------------------------------------#
cinch_load_extras(LEGION)

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
cinch_load_extras(MPI)
list (APPEND FleCSPH_LIBRARIES ${CINCH_RUNTIME_LIBRARIES} )

#------------------------------------------------------------------------------#
# Add OpenMP
#------------------------------------------------------------------------------#
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#------------------------------------------------------------------------------#
# Add GSL
#------------------------------------------------------------------------------#
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND FleCSPH_LIBRARIES ${GSL_LIBRARIES})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS}")

#------------------------------------------------------------------------------#
# Add Boost
#------------------------------------------------------------------------------#
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

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
# Add mpisph tests
#------------------------------------------------------------------------------#
cinch_add_application_directory("mpisph/")
cinch_add_application_directory("include/physics/test")
cinch_add_application_directory("include/tree_topology/test")

#------------------------------------------------------------------------------#
# Add application targets
#------------------------------------------------------------------------------#
cinch_add_application_directory("app/id_generators")
cinch_add_application_directory("app/drivers")
