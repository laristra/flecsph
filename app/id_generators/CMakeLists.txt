#----------------------------------------------------------------------------#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Executables
#------------------------------------------------------------------------------#

# see 'drivers/CMakeLists.txt` for explaination
# very similar function to driver, may consider
# generalizing
function(add_generator generator_name dim_list)
  foreach(dim ${dim_list})
    set(exe_name "${generator_name}_${dim}d_generator")

    add_executable(${exe_name})

    target_sources(${exe_name}
      PRIVATE
        ${generator_name}/main.cc
    )

    target_compile_definitions(${exe_name}
      PRIVATE
        EXT_GDIMENSION=${dim}
    )

    # NOTE: using imported 'FleCSI::flecsi' will
    # also append runtime_driver.cc, causing unresolved references
    # at linktime.
    # Use the explicit include to avoid this
    target_include_directories(${exe_name}
      PRIVATE
        ${CMAKE_SOURCE_DIR}/app/id_generators/include
        ${FleCSI_INCLUDE_DIR}
    )

    target_link_libraries(${exe_name}
      PRIVATE
        flecsph::flags
    )

    install(TARGETS ${exe_name} RUNTIME DESTINATION bin)
  endforeach()
endfunction()

add_generator(sodtube "1;2;3")
add_generator(sedov "1;2;3")
add_generator(noh "1;2;3")
add_generator(implosion "2;3")
add_generator(collapse "3")
add_generator(wtunnel "2;3")
add_generator(KH "2;3")
add_generator(RT "2;3")
add_generator(fluid "2;3")

# see 'drivers/CMakeLists.txt` for explaination
# very similar function to driver, may consider
# generalizing
function(add_generator_test generator_name partile_name parfile)
  set(TEST_NAME "${generator_name}_${partile_name}_test")
  #add_executable(${TEST_NAME} ALIAS ${generator_name})
  add_custom_target(${TEST_NAME} ${generator_name})

  configure_file(${PROJECT_SOURCE_DIR}/data/${parfile} "${CMAKE_BINARY_DIR}/tests" COPYONLY)
  add_test(
    NAME ${TEST_NAME}
    COMMAND ${MPIEXEC} -n 1 ${CMAKE_CURRENT_BINARY_DIR}/${generator_name} ${PROJECT_SOURCE_DIR}/data/${parfile}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
  )
  set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
  set_target_properties(${TEST_NAME} PROPERTIES FOLDER tests)
endfunction()

# add tests
add_generator_test(collapse_3d_generator collapse collapse_nx10.par)
add_generator_test(sodtube_1d_generator sodtube sodtube_t1_n100.par)
add_generator_test(RT_2d_generator RT RT_2d.par)
add_generator_test(KH_2d_generator KH KH_2d.par)
add_generator_test(sedov_2d_generator sedov sedov_nx20.par)
add_generator_test(sedov_3d_generator mesa mesa_nx20.par)
add_generator_test(sedov_3d_generator wvt wvt_nx20.par)
add_generator_test(noh_2d_generator noh noh_nx20.par)
add_generator_test(implosion_2d_generator implosion implosion_nx20.par)

