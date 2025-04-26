cmake_minimum_required(VERSION 2.8)

project(zjucad-ptree)

find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
  message("--[info] can not find Boost")
endif()

file(GLOB SOURCE *.cpp)
add_library(zjucad-ptree SHARED ${SOURCE} ${HEADER})
