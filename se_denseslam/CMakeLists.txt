cmake_minimum_required(VERSION 3.10)

project(se-denseslam)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

if(APPLE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__APPLE__")
endif(APPLE)

find_package(Eigen3 REQUIRED)
find_package(Sophus REQUIRED)
find_package(TooN REQUIRED)
if(WITH_OPENMP)
  find_package(OpenMP)
endif()

set(appname ${PROJECT_NAME})


# ---- PREPARE COMMON DEPENDENCIES  ------------ 

include_directories(./include ${TOON_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} 
  ${SOPHUS_INCLUDE_DIR} ../se_shared ./include
  ../se_core/include/ ../utils)
set(common_compile_flags "")
set(common_libraries lodepng)

# ----------------- OFUsion -----------------
set(field_type FIELD_TYPE=Occupancy)

add_library(${appname}-ofusion STATIC ../../se_denseslam/src/DenseSLAMSystem.cpp)
if(WITH_OPENMP AND OPENMP_FOUND)
  message(STATUS "Compiling with OpenMP support")
  target_compile_options(${appname}-ofusion PRIVATE ${OpenMP_CXX_FLAGS} -Wall -Wextra -Werror -Wno-unknown-pragmas)
  target_compile_definitions(${appname}-ofusion PRIVATE ${field_type})
  target_link_libraries(${appname}-ofusion   ${common_libraries})	
else()
  target_compile_options(${appname}-ofusion PRIVATE -Wall -Wextra -Werror -Wno-unknown-pragmas)
  target_compile_definitions(${appname}-ofusion PRIVATE ${field_type})
  target_link_libraries(${appname}-ofusion   ${common_libraries})	
endif()

list(APPEND BUILT_LIBS ${appname}-ofusion)

# ----------------- SDF -----------------
set(field_type FIELD_TYPE=SDF)

add_library(${appname}-sdf  ../../se_denseslam/src/DenseSLAMSystem.cpp)
if(WITH_OPENMP AND OPENMP_FOUND)
  message(STATUS "Compiling with OpenMP support")
  target_compile_options(${appname}-sdf PRIVATE ${OpenMP_CXX_FLAGS} -Wall -Wextra -Werror -Wno-unknown-pragmas)
  target_compile_definitions(${appname}-sdf PRIVATE ${field_type})
  target_link_libraries(${appname}-sdf ${OpenMP_CXX_FLAGS} ${common_libraries})	
else()
  target_compile_options(${appname}-sdf PRIVATE -Wall -Wextra -Werror -Wno-unknown-pragmas)
  target_compile_definitions(${appname}-sdf PRIVATE ${field_type})
  target_link_libraries(${appname}-sdf   ${common_libraries})	
endif()

list(APPEND BUILT_LIBS ${appname}-sdf)

set(BUILT_LIBS ${BUILT_LIBS} PARENT_SCOPE)