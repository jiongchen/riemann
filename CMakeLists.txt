project(riemann)
cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -fpermissive -fopenmp -w")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/external)
include_directories(${CMAKE_SOURCE_DIR}/external/eigen)

add_subdirectory(${CMAKE_SOURCE_DIR}/external/jtflib/util)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/jtflib/mesh)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/zjucad/ptree)

# BOOST
find_package(Boost COMPONENTS filesystem system program_options REQUIRED)
if(Boost_FOUND)
  message("-- Boost @ ${Boost_INCLUDE_DIRS}")
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

# LBFGS
find_package(LBFGS REQUIRED)
if(LBFGS_FOUND)
  message("-- LBFGS @ ${LBFGS_INCLUDE_DIRS}")
  include_directories(${LBFGS_INCLUDE_DIRS})
endif(LBFGS_FOUND)

# ALGLIB
find_package(Alglib REQUIRED)
if(ALGLIB_FOUND)
  message("-- Alglib @ ${ALGLIB_INCLUDES}")
  include_directories(${ALGLIB_INCLUDES})
endif(ALGLIB_FOUND)

# UMFPACK
find_package(UMFPACK REQUIRED)
if(UMFPACK_FOUND)
  include_directories(${UMFPACK_INCLUDES})
endif(UMFPACK_FOUND)

# SUITESPARSE
find_package(SuiteSparse REQUIRED)
if(SuiteSparse_FOUND)
  message("-- SuiteSparse @ ${SuiteSparse_INCLUDE_DIRS}")
  include_directories(${SuiteSparse_INCLUDE_DIRS})
endif(SuiteSparse_FOUND)

# MINPACK
find_package(MINPACK REQUIRED)
if(MINPACK_FOUND)
  message("-- minpack @ ${MINPACK_INCLUDE_DIRS}")
endif(MINPACK_FOUND)

# Arpack
find_package(Arpack REQUIRED)
if(ARPACK_FOUND)
  message("-- arpack @ ${ARPACK_LIBRARIES}")
endif(ARPACK_FOUND)

add_subdirectory(src)
add_subdirectory(examples)
