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
# Find Legion
#------------------------------------------------------------------------------#

if(FLECSI_RUNTIME_MODEL STREQUAL "legion" OR 
    FLECSI_RUNTIME_MODEL STREQUAL "mpilegion")

  find_package(Legion REQUIRED)

  message(STATUS "Legion found: ${Legion_FOUND}")

endif()

#
# Serial interface
#
if(FLECSI_RUNTIME_MODEL STREQUAL "serial")
  set(_runtime_path ${PROJECT_SOURCE_DIR}/flecsi/execution/serial)

#
# Legion interface
#
elseif(FLECSI_RUNTIME_MODEL STREQUAL "legion")
  set(_runtime_path ${PROJECT_SOURCE_DIR}/flecsi/execution/legion)
  if(NOT APPLE)
    set(FLECSI_RUNTIME_LIBRARIES  -ldl ${Legion_LIBRARIES} ${MPI_LIBRARIES})
  else()
    set(FLECSI_RUNTIME_LIBRARIES  ${Legion_LIBRARIES} ${MPI_LIBRARIES})
  endif()
    include_directories(${Legion_INCLUDE_DIRS})

#
# MPI interface
#
elseif(FLECSI_RUNTIME_MODEL STREQUAL "mpi")
  set(_runtime_path ${PROJECT_SOURCE_DIR}/flecsi/execution/mpi)
  if(NOT APPLE)
    set(FLECSI_RUNTIME_LIBRARIES  -ldl ${MPI_LIBRARIES})
  else()
    set(FLECSI_RUNTIME_LIBRARIES ${MPI_LIBRARIES})
  endif()

#
# MPI+Legion interface
#
elseif(FLECSI_RUNTIME_MODEL STREQUAL "mpilegion")
  if(NOT ENABLE_MPI)
    message (FATAL_ERROR "MPI is required for the mpilegion runtime model")
  endif()                    
    set(_runtime_path ${PROJECT_SOURCE_DIR}/flecsi/execution/mpilegion)
  if(NOT APPLE)
    set(FLECSI_RUNTIME_LIBRARIES  -ldl ${Legion_LIBRARIES} ${MPI_LIBRARIES})
  else()
    set(FLECSI_RUNTIME_LIBRARIES  ${Legion_LIBRARIES} ${MPI_LIBRARIES})
  endif()
  include_directories(${Legion_INCLUDE_DIRS})

#
# Default
#
else(FLECSI_RUNTIME_MODEL STREQUAL "serial")
  message(FATAL_ERROR "Unrecognized runtime selection")  
endif(FLECSI_RUNTIME_MODEL STREQUAL "serial")
