###############################################
# SetupFlags
#
# Purpose: top-level compiler/linker flags
###############################################

#------------------------------------------------
# more readable generator expressions
#------------------------------------------------
set(debug_tree "$<BOOL:${ENABLE_DEBUG_TREE}>")
set(build_debug "$<CONFIG:Debug>")
set(build_release "$<CONFIG:Release>")
set(unit_tests "$<BOOL:${ENABLE_UNIT_TESTS}>")
set(sys_cray "$<PLATFORM_ID:CrayLinuxEnvironment>")
set(cxx_intel "$<COMPILE_LANG_AND_ID:CXX,Intel>")
set(cxx_gnu "$<COMPILE_LANG_AND_ID:CXX,GNU>")
set(cxx_cray "$<COMPILE_LANG_AND_ID:CXX,Cray>")

# set C++17
target_compile_features(flecsph::flags
    INTERFACE
        cxx_std_17)

# some general definitions; subdirectories may define
# targets with specialized definitions
target_compile_definitions(flecsph::flags
    INTERFACE
        "LOG_STRIP_LEVEL=${LOG_STRIP_LEVEL}"
        "PARALLEL_IO"
        $<${debug_tree}:
          "ENABLE_DEBUG_TREE"
        >
)

# compiler-specific flags
# TODO: future, may be moved to it's own module
target_compile_options(flecsph::flags
    INTERFACE
      $<${build_debug}:
        "-g;-O2"
        $<${cxx_intel}:
          "-traceback"
        >
      >
      $<${build_release}:
        "-Ofast"
        $<${cxx_gnu}:
          "-ftree-vectorize"
        >
        # cray builds are cross-platorm, so don't
        # build with host default.
        # cray compiler wrapper _should_ enforce this
        # but just to make sure.
        $<$<AND:${cxx_intel},${sys_cray}>:
         "-xCORE-AVX2"
        >
        $<$<AND:${cxx_gnu},${sys_cray}>:
          "-march=core-avx2"
        >
      >
)

target_link_options(flecsph::flags
    INTERFACE
      # mpi issues with ipo
      $<${cxx_intel}:
        "-no-ipo"
      >
)

# global includes
target_include_directories(flecsph::flags
    INTERFACE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/include/physics
        ${CMAKE_SOURCE_DIR}/include/physics/eos
        ${CMAKE_SOURCE_DIR}/mpisph
        ${CMAKE_SOURCE_DIR}/app/drivers/include
)

# global libraries
# NOTE: imported libraries bring in includes, definitions, libs; convienent!
target_link_libraries(flecsph::flags
    INTERFACE
        Threads::Threads
        OpenMP::OpenMP_CXX
        MPI::MPI_CXX
        GSL::gsl
        Boost::headers
        m
        $<${unit_tests}:
          "GTest::GTest"
          "GTest::Main"
        >
        ${HDF5_LIBRARIES}
)

# HDF5 doesn't provide imported interface (as far as I can tell),
# so explicitily provide
target_include_directories(flecsph::flags
    INTERFACE
        ${HDF5_INCLUDE_DIR}
)
