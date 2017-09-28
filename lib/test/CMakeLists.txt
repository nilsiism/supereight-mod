cmake_minimum_required(VERSION 3.2)
project(octree_lib_unit_testing)

# GTest Root - Change to reflect your install dir
set(GTEST_ROOT ~/software/googletest/googletest)
find_package(GTest REQUIRED)

enable_testing()
add_compile_options(-std=c++1y )
include_directories(../../lib ../../utils)

add_subdirectory(allocation)
add_subdirectory(interp)
add_subdirectory(multiscale)
add_subdirectory(algorithms)
add_subdirectory(geometry)
