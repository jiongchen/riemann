# - Try to find PETSc
# Once done this will define
#
# PETSC_FOUND - system has PETSc
# PETSC_INCLUDES - the PETSc include directories
# PETSC_LIBRARIES - Link these to use PETSc
# PETSC_COMPILER - Compiler used by PETSc, helpful to find a compatible MPI
# PETSC_DEFINITIONS - Compiler switches for using PETSc
# PETSC_MPIEXEC - Executable for running MPI programs
# PETSC_VERSION - Version string (MAJOR.MINOR.SUBMINOR)
#
# Usage:
# find_package (PETSc) - find the directory a
#
# Setting these changes the behavior of the search
# PETSC_DIR - directory in which PETSc resides
# PETSC_ARCH - build architecture
#

# set (PETSC_DIR_ENV $ENV{PETSC_DIR})
# set (PETSC_ARCH_ENV $ENV{PETSC_ARCH})

# message (STATUS "Env. var. PETSC_DIR: $ENV{PETSC_DIR}")
# message (STATUS "Env. var. PETSC_ARCH: $ENV{PETSC_ARCH}")

find_path (PETSC_INCLUDE_DIR petsc.h
  HINTS ENV PETSC_DIR
  PATHS
  ${PETSC_DIR}/include
  # Debian paths
  /usr/lib/petscdir/3.4
  /usr/lib/petscdir/3.3 /usr/lib/petscdir/3.2 /usr/lib/petscdir/3.1
  /usr/lib/petscdir/3.0.0 /usr/lib/petscdir/2.3.3 /usr/lib/petscdir/2.3.2
  $ENV{HOME}/petsc
  DOC "PETSc Directory")

find_path (PETSC_INCLUDE_DIR_CONF petscconf.h
	   PATHS
	   ${PETSC_DIR}/${PETSC_ARCH}/include
	   DOC "PETSc Configuration Include")

find_library(PETSC_LIBRARIES 
             NAMES petsc
             PATHS
	      ${PETSC_DIR}/${PETSC_ARCH}/lib
	     DOC "PETSc libraries")

if (PETSC_LIBRARIES)
   set (PETSC_INCLUDE_DIR 	${PETSC_INCLUDE_DIR} ${PETSC_INCLUDE_DIR_CONF})
   set (PETSC_LIBRARY_DIR 	${PETSC_DIR}/${PETSC_ARCH}/lib)
   set (PETSC_FOUND ON)
endif (PETSC_LIBRARIES)


include(FindPackageHandleStandardArgs)

find_package_handle_standard_args (PETSc
  "PETSc could not be found. Be sure to set PETSC_DIR and PETSC_ARCH."
  PETSC_INCLUDE_DIR PETSC_LIBRARIES PETSC_LIBRARY_DIR)
