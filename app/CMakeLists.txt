#------------------------------------------------------------------------------#
# Debug and release flags     
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG "-O2 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wno-unused -Wno-comment") 

# includes 
include_directories(${CMAKE_SOURCE_DIR}/include)

#------------------------------------------------------------------------------#
# Add REQUIRED PACKAGES
#------------------------------------------------------------------------------#

# Pthreads for parallelism in thre tree 
find_package(Threads REQUIRED)
include_directories(${THREADS_INCLUDE_PATH})

# Add MPI and Legion
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(Legion REQUIRED)
include_directories(${Legion_INCLUDE_DIRS})

# Add FleCSI TODO= Change findFlecsi.cmake
find_package(FleCSI REQUIRED)
include_directories(${FleCSI_INCLUDE_DIRS})

# Add OpenMP 
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#------------------------------------------------------------------------------#
# Add Not REQUIRED PACKAGES
#------------------------------------------------------------------------------#

# Add HDF5
find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIR})

#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# MPILEGION version
add_executable(tree_mpilegion 
  tree.cc 
  tree_driver.cc 
  io.cc 
  physics.cc
  /home/jloiseau/flecsi/flecsi/flecsi/execution/mpilegion/runtime_driver.cc
)

target_link_libraries(tree_mpilegion ${FleCSI_LIBRARIES}) 
target_link_libraries(tree_mpilegion ${LEGION_LIBRARY} ${REALM_LIBRARY})
target_link_libraries(tree_mpilegion /home/jloiseau/local/lib64/liblegion.so rt dl)
target_link_libraries(tree_mpilegion ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(tree_mpilegion ${MPI_LIBRARIES})
target_link_libraries(tree_mpilegion ${HDF5_LIBRARIES})
#target_compile_definitions(tree_mpilegion PRIVATE -DFLECSI_RUNTIME_MODEL=FLECSI_RUNTIME_MODEL_mpi)
target_compile_definitions(tree_mpilegion PRIVATE -DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)
