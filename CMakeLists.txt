cmake_minimum_required(VERSION 3.10)

project(se)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

add_compile_options(-std=c++14 -march=native)
add_subdirectory(se_tools)

# Build se_denseslam lib
find_package(Eigen3 REQUIRED)
find_package(Sophus REQUIRED)

option(WITH_OPENMP "Compile with OpenMP" ON)


set(BUILT_LIBS "")
add_subdirectory(se_denseslam)
add_subdirectory(se_apps)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/se_core/include/ DESTINATION include)
install(FILES ${CMAKE_SOURCE_DIR}/se_denseslam/include/se/volume_traits.hpp DESTINATION include/se)