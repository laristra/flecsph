#------------------------------------------------------------------------------#
# Find FleCSI, This should be in FindFleCSI
#------------------------------------------------------------------------------#
FIND_PATH(FleCSI_INCLUDE_DIR NAMES flecsi.h flecsi-config.h)
# Look for the library.
FIND_LIBRARY(FleCSI_LIBRARY NAMES flecsi libflecsi FleCSI libFleCSI)

# Look for the runtime driver
FIND_PATH(FleCSI_RUNTIME
	NAMES
		runtime_driver.cc
	PATH_SUFFIXES
		share/FleCSI/runtime share/flecsi/runtime
)

# handle the QUIETLY and REQUIRED arguments and set FleCSI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FleCSI DEFAULT_MSG FleCSI_INCLUDE_DIR FleCSI_LIBRARY FleCSI_RUNTIME)

if(FleCSI_FOUND)
    # Copy the results to the output variables.
    SET(FleCSI_INCLUDE_DIRS ${FleCSI_INCLUDE_DIR})
    SET(FleCSI_LIBRARIES ${FleCSI_LIBRARY})
    SET(FleCSI_RUNTIME ${FleCSI_RUNTIME})
else()
    SET(FleCSI_INCLUDE_DIRS)
    SET(FleCSI_LIBRARIES)
    SET(FleCSI_RUNTIME)
endif()

add_library(FleCSI::flecsi UNKNOWN IMPORTED)
set_target_properties(FleCSI::flecsi PROPERTIES
    IMPORTED_LOCATION "${FleCSI_LIBRARY}"
    INTERFACE_SOURCES "${FleCSI_RUNTIME}/runtime_driver.cc"
    INTERFACE_INCLUDE_DIRECTORIES "${FleCSI_INCLUDE_DIR}"
)

MARK_AS_ADVANCED(FleCSI_INCLUDE_DIR FleCSI_LIBRARY FleCSI_RUNTIME)

