cmake_minimum_required(VERSION 2.8)
project(jtf-util)

file(GLOB source *.cpp)

add_library(jtf-util SHARED ${header} ${source})
