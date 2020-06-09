#------------------------------------------------------------------------------#
# Copyright (c) 2020 Triad National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# The project name
#------------------------------------------------------------------------------#
project(flecsph)
cmake_minimum_required(VERSION 3.15)

#------------------------------------------------
# kick the tires
#------------------------------------------------

enable_language(CXX)
include(CMakeDetermineCXXCompiler)
include(CheckCXXCompilerFlag)

#------------------------------------------------
# global interface targets
#------------------------------------------------

# interface libraries
# see: cmake/SetupDeps.cmake, cmake/SetupFlags.cmake
#add_library(flecsph::include_sources INTERFACE IMPORTED)
#add_library(flecsph::compile_flags INTERFACE IMPORTED)
#add_library(flecsph::library_flags INTERFACE IMPORTED)
add_library(flecsph::flags INTERFACE IMPORTED)

#------------------------------------------------
# construct configuration state
#------------------------------------------------

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# make sure order is correct when adding modules
include(ProcessOptions)
include(SetupDeps)
include(SetupFlags)
include(SetupTesting)

#add_subdirectory("include")
#------------------------------------------------------------------------------#
# Add mpisph tests
#------------------------------------------------------------------------------#
add_subdirectory("include")
add_subdirectory("mpisph/")
add_subdirectory("include/physics/test")
add_subdirectory("include/tree_topology/test")

#------------------------------------------------------------------------------#
# Add application targets
#------------------------------------------------------------------------------#
add_subdirectory("app/id_generators")
add_subdirectory("app/drivers")

#------------------------------------------------------------------------------#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#------------------------------------------------------------------------------#
