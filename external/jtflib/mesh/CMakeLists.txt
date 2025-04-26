cmake_minimum_required(VERSION 2.8)

project(jtf-mesh)

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

file(GLOB source *.cpp)

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
    message("CAN NOT FIND BOOST")
endif(Boost_FOUND)

add_library(jtf-mesh SHARED ${source})
target_link_libraries(jtf-mesh
  jtf-util)
