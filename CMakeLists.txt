project(riemann)
cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -fpermissive -fopenmp")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")
include(geo_sim_sdk)
include_geo_sim_sdk()
link_geo_sim_sdk()
include_directories(${CMAKE_SOURCE_DIR})

# BOOST
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS filesystem system REQUIRED)

# EIGEN3
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

# UMFPACK
find_package(UMFPACK REQUIRED)
if(UMFPACK_FOUND)
    include_directories(${UMFPACK_INCLUDES})
endif(UMFPACK_FOUND)

# SUITESPARSE
find_package(SuiteSparse REQUIRED)
if(SuiteSparse_FOUND)
    include_directories(${SuiteSparse_INCLUDE_DIRS})
endif(SuiteSparse_FOUND)

# IPOPT
find_package(IPOPT REQUIRED)
if(IPOPT_FOUND)
    message("IPOPT found")
    include_directories(${IPOPT_INCLUDE_DIRS})
endif(IPOPT_FOUND)

# ViennaCL
find_package(ViennaCL)
if(VIENNACL_FOUND)
    message("ViennaCL found")
    include_directories(${VIENNACL_INCLUDE_DIRS})
endif(VIENNACL_FOUND)

add_subdirectory(src)
add_subdirectory(examples)
