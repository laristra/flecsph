CMAKE_MINIMUM_REQUIRED(VERSION 3.3)

PROJECT(ScalingFramework CXX)

# Compiler settings
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS "-fPIC")
set(CMAKE_CXX_STANDARD 11)


# Find Boost
find_package(Boost COMPONENTS python REQUIRED)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIR})
elseif (NOT Boost_FOUND)
  	message(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
endif()


# Find Python
find_package( PythonLibs REQUIRED )
if (PythonLibs_FOUND)
	include_directories(${PYTHON_INCLUDE_DIR})
	message(STATUS "${PYTHON_INCLUDE_DIR}")
elseif (NOT PythonLibs 2.7)
	message(FATAL_ERROR "Python not found")
endif()



# Include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/IOTests)
include_directories(${CMAKE_SOURCE_DIR}/utils)
include_directories(/home/pascal/NGC/IO/phdf5-1.8.18/install/include/)
include_directories(/home/pascal/NGC/IO/hdf5-1.8.18/install/include/)
include_directories(/home/pascal/NGC/IO/H5hut-1.99.13/install/include/)

# Link directories
link_directories(${CMAKE_SOURCE_DIR}/build)
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


# Library
add_library(hdf5ParticleIO SHARED IOTests/hdf5ParticleIO.h)
set_target_properties(hdf5ParticleIO PROPERTIES LINKER_LANGUAGE "CXX")
target_link_libraries(hdf5ParticleIO hdf5.so H5hut.a)
#set_property(TARGET hdf5ParticleIO PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES hdf5.a H5hut.a)

# Python 
add_library(hdf5PythonWrapper SHARED IOTests/hdf5PythonWrapper.cpp)
target_link_libraries( hdf5PythonWrapper ${Boost_LIBRARIES} hdf5ParticleIO.so hdf5.so H5hut.so)
set_target_properties( hdf5PythonWrapper PROPERTIES PREFIX "" )


# Executable
add_executable(ScalingFramework main.cpp)
target_link_libraries(ScalingFramework  hdf5.so H5hut.a)
