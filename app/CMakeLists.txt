#------------------------------------------------------------------------------#
# Debug and release flags     
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG 
  "-O2 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")
set(CMAKE_CXX_FLAGS_RELEASE 
  "-O3 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses") 

# includes 
include_directories(${CMAKE_SOURCE_DIR}/include)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# MPILEGION version
add_executable(tree_mpilegion 
  main.cc 
  main_driver.cc 
  mpi_partition.cc
  io.cc
  sodtube.cc 
  physics.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(tree_mpilegion ${FleCSPH_LIBRARIES})

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_DRIVER)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)
