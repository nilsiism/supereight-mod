cmake_minimum_required(VERSION 2.8)
project(octree_lib_unit_testing)

# GTest Root - Change to reflect your install dir
set(GTEST_ROOT ~/software/googletest/googletest)
find_package(GTest REQUIRED)

enable_testing()
include_directories(../../lib ../../utils)

add_subdirectory(allocation)
