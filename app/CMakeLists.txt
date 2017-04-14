#------------------------------------------------------------------------------#
# Debug and release flags     
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG "-O2 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall") 

# includes 
include_directories(${CMAKE_SOURCE_DIR}/include)

#------------------------------------------------------------------------------#
# Threads needed for this version
#------------------------------------------------------------------------------#

# Pthreads for parallelism in thre tree 
find_package(Threads REQUIRED)
include_directories(${THREADS_INCLUDE_PATH})

# Add MPI and Legion
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(Legion REQUIRED)
include_directories(${Legion_INCLUDE_DIRS})

find_package(FleCSI REQUIRED)
include_directories(${FleCSI_INCLUDE_DIRS})

message(STATUS "Includes: ${FleCSI_INCLUDE_DIRS}")
message(STATUS "Libraries ${FleCSI_LIBRARIES}")

#------------------------------------------------------------------------------#
# Add flecsi with absolute path and link library    
#------------------------------------------------------------------------------#

# For flecsi.h
#include_directories(/home/jloiseau/local/include)

#include_directories(/home/jloiseau/local/include/legion)
#include_directories(/home/jloiseau/local/include/realm)
#include_directories(/home/jloiseau/local/include/mappers)

#include_directories(/home/jloiseau/local/include/flecsi)

# for legion/helper.h 
#include_directories(/home/jloiseau/flecsi/flecsi/)


#------------------------------------------------------------------------------#
# Executables    
#------------------------------------------------------------------------------#

# Serial version
#add_executable(tree_serial tree.cc tree_driver.cc io.cc physics.cc)
#target_link_libraries(tree_serial /home/jloiseau/flecsi/flecsi/build/lib/libflecsi.so)
#target_link_libraries(tree_serial ${CMAKE_THREAD_LIBS_INIT})

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
#target_compile_definitions(tree_mpilegion PRIVATE -DFLECSI_RUNTIME_MODEL=FLECSI_RUNTIME_MODEL_mpi)
target_compile_definitions(tree_mpilegion PRIVATE -DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)
