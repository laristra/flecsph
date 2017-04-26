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
include_directories(${CMAKE_SOURCE_DIR}/utils)


# Sources
set(SRC_LIST main.cpp)


# Options
option (ENABLE_DEBUG "enable debugging mode" ON)


# Executable
add_executable(ScalingFramework main.cpp)
