#------------------------------------------------------------------------------#
# Debug and release flags     
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG 
  "-g -O2 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")
set(CMAKE_CXX_FLAGS_RELEASE 
  "-O3 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses") 

# includes 
include_directories(${CMAKE_SOURCE_DIR}/include)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# MPILEGION version
add_library(mpisph 
  mpi_partition.cc
  io.cc
  #${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(mpisph ${FleCSPH_LIBRARIES})

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_DRIVER)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)
