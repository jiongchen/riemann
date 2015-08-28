#
# Try to find ViennaCL library and include path.
# Once done this will define
#
# VIENNACL_FOUND
# VIENNACL_INCLUDE_DIRS
# VIENNACL_LIBRARIES
# 

option(ENABLE_CUDA "Use the CUDA backend" ON)
option(ENABLE_OPENMP "Use OpenMP acceleration" ON)
option(ENABLE_EIGEN "Use Eigen" ON)
option(ENABLE_OPENCL "Use the OpenCL backend" OFF)
option(ENABLE_ASAN "Build with address sanitizer if available" OFF)

if(ENABLE_EIGEN)
    message("ViennaCL enables Eigen")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVIENNACL_WITH_EIGEN -DVIENNACL_HAVE_EIGEN")
endif(ENABLE_EIGEN)

if(ENABLE_OPENCL)
  find_package(OpenCL REQUIRED)
endif(ENABLE_OPENCL)

if(ENABLE_CUDA)
   find_package(CUDA REQUIRED)
   set(CUDA_ARCH_FLAG "-arch=sm_20" CACHE STRING "Use one out of sm_13, sm_20, sm_30, ...")
   set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "${CUDA_ARCH_FLAG}" "-DVIENNACL_WITH_CUDA" "-std=c++11")
endif(ENABLE_CUDA)

if(ENABLE_OPENMP)
   find_package(OpenMP REQUIRED)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DVIENNACL_WITH_OPENMP")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DVIENNACL_WITH_OPENMP")
   set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    ${OpenMP_EXE_LINKER_FLAGS}")
   set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${OpenMP_MODULE_LINKER_FLAGS}")
   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
   set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${OpenMP_STATIC_LINKER_FLAGS}")
endif(ENABLE_OPENMP)

if(ENABLE_ASAN)
  add_c_compiler_flag_if_supported("-fsanitize=address")
  add_c_linker_flag_if_supported("-fsanitize=address")
endif(ENABLE_ASAN)

find_path(VIENNACL_INCLUDE_DIR viennacl/forwards.h
  PATHS /usr/include
  DOC "The ViennaCL include path"
)

include(FindPackageHandleStandardArgs)
if(ENABLE_CUDA)
  set(VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
  set(VIENNACL_LIBRARIES ${CUDA_LIBRARIES})
  find_package_handle_standard_args(ViennaCL "ViennaCL not found!" VIENNACL_INCLUDE_DIR CUDA_INCLUDE_DIRS CUDA_LIBRARIES)
else(ENABLE_CUDA)
  set(VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIR})
  set(VIENNACL_LIBRARIES "")
  find_package_handle_standard_args(ViennaCL "ViennaCL not found!" VIENNACL_INCLUDE_DIR)
endif(ENABLE_CUDA)

mark_as_advanced(VIENNACL_INCLUDE_DIR)
