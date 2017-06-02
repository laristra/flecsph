CMAKE_MINIMUM_REQUIRED(VERSION 3.3)

PROJECT(ScalingFramework CXX)

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
include_directories(/home/pascal/software/phdf5-1.8.18/install/include/)
include_directories(/home/pascal/software/hdf5-1.8.18/install/include/)
include_directories(/home/pascal/software/H5hut-1.99.13/install/include/)

# Link directories
link_directories(/home/pascal/software/phdf5-1.8.18/install/lib/)
link_directories(/home/pascal/software/hdf5-1.8.18/install/lib/)
link_directories(/home/pascal/software/H5hut-1.99.13/install/lib/)


# Sources
set(SRC_LIST main.cpp)


# Options
option (ENABLE_DEBUG "enable debugging mode" ON)
if (ENABLE_DEBUG)
	add_definitions(-DENABLE_DEBUG=1)
endif()


#Library
add_library(h5FlecsiPython SHARED IOTests/hdf5SimIO.h)
set_target_properties(h5FlecsiPython PROPERTIES LINKER_LANGUAGE "CXX")
set_property(TARGET h5FlecsiPython PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES hdf5.so H5hut.a)


# Executable
add_executable(ScalingFramework main.cpp)

target_link_libraries(ScalingFramework hdf5.so H5hut.a)
