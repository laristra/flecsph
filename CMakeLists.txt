CMAKE_MINIMUM_REQUIRED(VERSION 3.3)

PROJECT(ScalingFramework)

# Compiler settings
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS "-fPIC")
set(CMAKE_CXX_STANDARD 11)


# Include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/IOTests)
include_directories(${CMAKE_SOURCE_DIR}/utils)
include_directories(/home/pascal/NGC/IO/phdf5-1.8.18/install/include/)
include_directories(/home/pascal/NGC/IO/hdf5-1.8.18/install/include/)
include_directories(/home/pascal/NGC/IO/H5hut-1.99.13/install/include/)

# Link directories
link_directories(/home/pascal/NGC/IO/phdf5-1.8.18/install/lib/)
link_directories(/home/pascal/NGC/IO/hdf5-1.8.18/install/lib/)
link_directories(/home/pascal/NGC/IO/H5hut-1.99.13/install/lib/)


# Sources
set(SRC_LIST main.cpp)


# Options
option (ENABLE_DEBUG "enable debugging mode" ON)
if (ENABLE_DEBUG)
	add_definitions(-DENABLE_DEBUG=1)
endif()


# Executable
add_executable(ScalingFramework main.cpp)
target_link_libraries(ScalingFramework hdf5.so hdf5_cpp.so H5hut.a)
