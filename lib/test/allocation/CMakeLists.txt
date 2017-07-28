cmake_minimum_required(VERSION 2.8)
project(octree_lib)

enable_testing()
set(GTEST_ROOT ~/software/googletest/googletest)
find_package(GTest REQUIRED)

add_compile_options(-std=c++1y )

if(GTEST_FOUND)
  message("GTEST FOUND")
endif()

set(PROJECT_TEST_NAME alloc_unittest)
add_executable(${PROJECT_TEST_NAME} alloc_unittest.cpp)
target_include_directories(${PROJECT_TEST_NAME} PUBLIC ${GTEST_INCLUDE_DIRS})
target_link_libraries(${PROJECT_TEST_NAME} ${GTEST_BOTH_LIBRARIES} pthread)

 GTEST_ADD_TESTS(alloc_unittest "" alloc_unittest.cpp)
