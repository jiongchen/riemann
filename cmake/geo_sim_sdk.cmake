set(CMAKE_DEBUG_POSTFIX "d")

# set OS, COMPILER BITS

set(OS ${CMAKE_SYSTEM_NAME})

if(NOT DEFINED BITS)
  if(${CMAKE_SIZEOF_VOID_P} MATCHES "8")
	set(BITS "64")
  else(${CMAKE_SIZEOF_VOID_P} MATCHES "8")
	set(BITS "32")
  endif(${CMAKE_SIZEOF_VOID_P} MATCHES "8")
endif(NOT DEFINED BITS)

if(MSVC)
  if(${MSVC_VERSION} MATCHES "1500")
	set(COMPILER "vc9sp1")
  endif(${MSVC_VERSION} MATCHES "1500")
  if(${MSVC_VERSION} MATCHES "1600")
	set(COMPILER "vc2010")
  endif(${MSVC_VERSION} MATCHES "1600")
else(MSVC)
  exec_program(${CMAKE_C_COMPILER} ARGS --version OUTPUT_VARIABLE CMAKE_C_COMPILER_VERSION)
  if(CMAKE_C_COMPILER_VERSION MATCHES ".*3\\.[0-9].*")
	set(COMPILER "gcc3")
  elseif(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9].*")
	set(COMPILER "gcc4")
  endif(CMAKE_C_COMPILER_VERSION MATCHES ".*3\\.[0-9].*")
endif(MSVC)

macro(include_geo_sim_sdk)
  include_directories(
	$ENV{HOME}/usr/include
	$ENV{HOME}/usr/${OS}/include
	$ENV{HOME}/usr/${OS}/${BITS}/include
	$ENV{HOME}/usr/${OS}/${BITS}/${COMPILER}/include  # indeed, its better to be ${COMPILER}/${COMPILER_VERSION}
	)
endmacro(include_geo_sim_sdk)

macro(link_geo_sim_sdk)
  link_directories(
	$ENV{HOME}/usr/${OS}/lib
	$ENV{HOME}/usr/${OS}/${BITS}/lib
	$ENV{HOME}/usr/${OS}/${BITS}/${COMPILER}/lib
	)
endmacro(link_geo_sim_sdk)

message(${OS}/${BITS}/${COMPILER})

macro(geo_sim_sdk_install_cxx)
  install(${ARGV}
	RUNTIME DESTINATION ${OS}/${BITS}/${COMPILER}/bin
	ARCHIVE DESTINATION ${OS}/${BITS}/${COMPILER}/lib
	LIBRARY DESTINATION ${OS}/${BITS}/${COMPILER}/lib
	)
endmacro(geo_sim_sdk_install_cxx)

macro(geo_sim_sdk_install_c)
  install(${ARGV}
	RUNTIME DESTINATION ${OS}/${BITS}/bin
	ARCHIVE DESTINATION ${OS}/${BITS}/lib
	LIBRARY DESTINATION ${OS}/${BITS}/lib
	)
endmacro(geo_sim_sdk_install_c)

macro(geo_sim_sdk_install_header)
install(DIRECTORY include/
  DESTINATION include/${ARGV}
  PATTERN .svn EXCLUDE
	PATTERN *~ EXCLUDE
	PERMISSIONS GROUP_READ WORLD_READ)
endmacro(geo_sim_sdk_install_header)

# some helper code
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

macro(aux_source_directory_with_headers dir src-variable hdr-variable)
  aux_source_directory(${dir} ${src-variable})
  file(GLOB ${hdr-variable} ${dir} *.h)
  source_group("Header Files" FILES ${hdr-variable})
endmacro(aux_source_directory_with_headers)

macro(add_library_from_dir name mode dir)
  aux_source_directory_with_headers(${dir} sources headers)
  add_library(${name} ${mode} ${sources} ${headers})
endmacro(add_library_from_dir name mode dir)

macro(add_executable_from_dir name dir)
  aux_source_directory_with_headers(${dir} sources headers)
  add_executable(${name} ${sources} ${headers})
endmacro(add_executable_from_dir)

if(WIN32)
  add_definitions(/DGLUT_DISABLE_ATEXIT_HACK)
  if(MSVC)
	# disable iterator checking and dll interface warning.
	add_definitions("/D_SCL_SECURE_NO_WARNINGS /wd4251")
  endif(MSVC)
endif(WIN32)

CONFIGURE_FILE(
  "$ENV{HOME}/usr/share/cmake/Modules/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

macro(add_uninstall_target)
  ADD_CUSTOM_TARGET(UNINSTALL
	"${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
endmacro(add_uninstall_target)
