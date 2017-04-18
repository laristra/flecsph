#------------------------------------------------------------------------------#
# Debug and release flags     
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG "-O2 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wno-unused -Wno-comment") 

# includes 
include_directories(${CMAKE_SOURCE_DIR}/include)

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# MPILEGION version
add_executable(tree_mpilegion 
  tree.cc 
  tree_driver.cc 
  io.cc 
  physics.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

target_link_libraries(tree_mpilegion ${FleCSPH_LIBRARIES})

target_compile_definitions(tree_mpilegion PRIVATE
  -DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)
