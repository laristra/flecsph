###############################################
# ProcessOptions
#
# Purpose: set up project build options
###############################################

# for dependent option macro
include(CMakeDependentOption)

#------------------------------------------------
# determine/set user/configuration options
#------------------------------------------------

# sets up building unit tests
option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
# enables debug messages from tree
option(ENABLE_DEBUG_TREE "Enable debug tree" OFF)
# TODO: get rid of this
#option(ENABLE_DEBUG "Compile in DEBUG mode" OFF)
# sets integrated log level (0 - none, X - ?)
set(LOG_STRIP_LEVEL 0 CACHE STRING "LOG Strip level")
# integer width for keys
set(KEY_INTEGER_TYPE "uint64_t" CACHE STRING "Type of integer used to generate keys")
set_property(CACHE KEY_INTEGER_TYPE PROPERTY STRINGS "uint32_t" "uint64_t" "uint128_t")

# tentative; color output at building
option(ENABLE_FORCE_COMPILE_COLORED "Forces build to use colorized output" ON)

# TODO: Release or Debug default?
# Set a default build type
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------
# prettier output options
#------------------------------------------------

#------------------------------------------------
# determine/set system information
#------------------------------------------------

# get machine information
cmake_host_system_information(RESULT FSPH_hostname QUERY FQDN)

# detect if in CrayPE
# this will default to 1 if CMake think's this is on Cray,
# but the user may override it. If not on Cray, the user
# is not presented with the options to edit
cmake_dependent_option(FSPH_USE_CRAY_LINUX
  "Use Cray environment in configuration"
  ON
  "DEFINED ENV{CRAYPE_VERSION}" OFF
)

#------------------------------------------------
# post operations
#------------------------------------------------

# create directory if tests
if(ENABLE_UNIT_TESTS)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
endif()

# When building on Cray, tell CMake to use
# compute-node configuration to build.
if(FSPH_USE_CRAY_LINUX)
  set(CMAKE_SYSTEM_NAME "CrayLinuxEnvironment")
endif()

if(ENABLE_FORCE_COMPILE_COLORED)
  add_compile_options(
    "$<$<CXX_COMPILER_ID:GNU>:-fdiagnostics-color=always>"
  )
endif()
