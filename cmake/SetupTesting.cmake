# Macro to add googletest
macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN} ${FleCSI_RUNTIME}/runtime_driver.cc)

    target_link_libraries(${TESTNAME}
        PRIVATE
            flecsph::flags
            FleCSI::flecsi
    )
    add_test(
        NAME ${TESTNAME}
        COMMAND ${TESTNAME}
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests"
    )
    #set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    set_target_properties(${TESTNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests")

endmacro()

# Macro to add googletest for MPI
macro(package_add_test_MPI TESTNAME)
    add_executable(${TESTNAME} ${ARGN} ${FleCSI_RUNTIME}/runtime_driver.cc)
    target_link_libraries(${TESTNAME}
        PRIVATE
            flecsph::flags
            FleCSI::flecsi
    )
    add_test(
        NAME ${TESTNAME}
        COMMAND ${MPIEXEC} -n 4 "${PROJECT_BINARY_DIR}/tests/${TESTNAME}"
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/tests"
    )
    #set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    set_target_properties(${TESTNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests")

endmacro()
