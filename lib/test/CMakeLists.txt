cmake_minimum_required(VERSION 3.10)
project(octree_lib_unit_testing)

# GTest Root - Change to reflect your install dir
set(GTEST_ROOT ~/software/googletest/googletest)
find_package(GTest REQUIRED)
find_package(Eigen3 REQUIRED)

enable_testing()
add_compile_options(-std=c++1y -O3)
include_directories(../../lib ../../utils ../../lib/geometry ${EIGEN3_INCLUDE_DIR})

add_subdirectory(octree)
add_subdirectory(allocation)
add_subdirectory(functor)
add_subdirectory(interp)
add_subdirectory(multiscale)
add_subdirectory(algorithms)
add_subdirectory(geometry)
add_subdirectory(utils)
