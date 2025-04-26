# FindMINPACK.cmake
# Find the MINPACK library
#
# This module defines:
#  MINPACK_FOUND        - True if MINPACK is found
#  MINPACK_INCLUDE_DIRS - The MINPACK include directories
#  MINPACK_LIBRARIES    - The MINPACK libraries
#
# You can specify custom MINPACK library locations using:
#  MINPACK_ROOT         - Root directory for MINPACK installation

# Try to find MINPACK in standard locations or user-specified location
find_path(MINPACK_INCLUDE_DIR
  NAMES minpack.h
  PATHS
    ${MINPACK_ROOT}/include
    /usr/include
    /usr/local/include
    /opt/local/include
  DOC "MINPACK include directory"
)

find_library(MINPACK_LIBRARY
  NAMES minpack
  PATHS
    ${MINPACK_ROOT}/lib
    /usr/lib
    /usr/local/lib
    /opt/local/lib
  DOC "MINPACK library"
)

# Set the MINPACK_FOUND variable
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MINPACK
  DEFAULT_MSG
  MINPACK_LIBRARY
  MINPACK_INCLUDE_DIR
)

# Set the output variables
if(MINPACK_FOUND)
  set(MINPACK_LIBRARIES ${MINPACK_LIBRARY})
  set(MINPACK_INCLUDE_DIRS ${MINPACK_INCLUDE_DIR})
endif()

# Hide internal variables from CMake GUI
mark_as_advanced(MINPACK_INCLUDE_DIR MINPACK_LIBRARY)
