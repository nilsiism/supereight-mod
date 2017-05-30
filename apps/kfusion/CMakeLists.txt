cmake_minimum_required(VERSION 2.8)

project(pamela_testcases)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ")
  add_definitions(-DDEBUG=1)
else(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3 ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 ")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

add_compile_options(-std=c++11)

if(APPLE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__APPLE__")
endif(APPLE)


# ---- PREPARE COMMON DEPENDENCIES  ------------ 


function(add_version appname version compile_flags link_flags)

  if(APPLE)
    set(main_common_libraries     ${link_flags} stdc++ pthread)
  else(APPLE)
    set(main_common_libraries     ${link_flags} rt stdc++ pthread)
  endif(APPLE)

  set(main_common_compile_flags ${compile_flags})

  if (OPENNI2_FOUND)
    include_directories(${OPENNI2_INCLUDE_PATHS})
    set(main_common_libraries ${main_common_libraries} ${OPENNI2_LIBRARIES})
    set(main_common_compile_flags "-DDO_OPENNI ${main_common_compile_flags}")
  endif(OPENNI2_FOUND)



  #  ------ BENCHMARK --------
  add_executable(${appname}-benchmark-${version} src/benchmark.cpp src/PowerMonitor.cpp)
  target_link_libraries(${appname}-benchmark-${version} ${appname}-${version} ${main_common_libraries})
  SET_TARGET_PROPERTIES(${appname}-benchmark-${version}
    PROPERTIES COMPILE_FLAGS "${main_common_compile_flags}")


  #  ------ MAIN --------
  if (GLUT_FOUND)
    add_executable(${appname}-main-${version} src/mainQt.cpp src/PowerMonitor.cpp src/reader.cpp)

    SET_TARGET_PROPERTIES(${appname}-main-${version}
      PROPERTIES COMPILE_FLAGS "${main_common_compile_flags}")


    if(APPLE)
      target_link_libraries(${appname}-main-${version} ${appname}-${version} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES}  ${main_common_libraries} )
    else(APPLE)
      target_link_libraries(${appname}-main-${version} ${appname}-${version} ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES}  ${main_common_libraries} GL )
    endif(APPLE)
  endif(GLUT_FOUND)
  #  ------ QT --------
  if (Qt5OpenGL_FOUND)

    file(GLOB qt_srcs qt/*)

    QT5_ADD_RESOURCES( ResourceSources  "qt/images.qrc")
    set(qt_srcs ${qt_srcs} ${ResourceSources} include/PowerMonitor.h ../../commons/perfstats.h)
    SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)

    include_directories(qt)
    add_definitions( -D_FILE_OFFSET_BITS=64 -fPIC )

    add_executable(${appname}-qt-${version} src/mainQt.cpp src/PowerMonitor.cpp ${qt_srcs} src/reader.cpp)

    ### qt5_use_modules is not available with cmake 2.8.7 ###
    #qt5_use_modules(qtgui-${version} Core Gui Widgets OpenGL PrintSupport)

    add_definitions(${Qt5Widgets_DEFINITIONS})
    add_definitions(${Qt5OpenGL_DEFINITIONS})
    add_definitions(${Qt5PrintSuppoer_DEFINITIONS})
    include_directories(${Qt5Widgets_INCLUDE_DIRS})
    include_directories(${Qt5OpenGL_INCLUDE_DIRS})
    include_directories(${Qt5PrintSupport_INCLUDE_DIRS})
    include_directories(${Qt5Gui_INCLUDE_DIRS})

    SET_TARGET_PROPERTIES(${appname}-qt-${version} PROPERTIES COMPILE_FLAGS "${main_common_compile_flags} -std=c++0x  -D__QT__")

    if(APPLE)
      target_link_libraries(${appname}-qt-${version} Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL Qt5::PrintSupport 
        ${main_common_libraries} ${appname}-${version}  ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES}   )  
    else(APPLE) 
      target_link_libraries(${appname}-qt-${version} Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL Qt5::PrintSupport 
        ${main_common_libraries} ${appname}-${version}  ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES}  GL )  
    endif(APPLE)


  endif(Qt5OpenGL_FOUND)

endfunction()

set(appname kfusion)

find_package(TooN   REQUIRED)
find_package(OpenCL )
find_package(CUDA   )
find_package(GLUT   )
find_package(OpenGL )
find_package(OpenNI2)
find_package(PAPI)

find_package(PkgConfig)
find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt5Gui)
find_package(Qt5OpenGL)
find_package(Qt5PrintSupport)
include_directories(./include ${TOON_INCLUDE_PATHS} ./thirdparty ../../lib ../../commons ${PCL_INCLUDE_DIRS})

# ---- PREPARE COMMON DEPENDENCIES  ------------ 
add_library(tiny_obj thirdparty/tiny_obj_loader.cc)
if(APPLE)
  set(common_libraries stdc++ pthread tiny_obj lodepng)
else(APPLE)
  set(common_libraries rt stdc++ pthread tiny_obj lodepng)
endif(APPLE)


set(common_compile_flags "")

OPTION(STATS "Enable stats collection" OFF)
IF(STATS)
  add_definitions(-DSTATS)
ENDIF(STATS)

# ----------------- CPP VERSION ----------------- 

add_library(${appname}-cpp  src/cpp/kernels.cpp)
target_link_libraries(${appname}-cpp   ${common_libraries})	
add_version(${appname} cpp "" "")

# # ----------------- CPP VERSION OCTREE ----------------- 
# 
# add_library(${appname}-octree-cpp  src/octree/kernels.cpp src/octree/data-structure/octree.hpp)
# target_link_libraries(${appname}-octree-cpp  ${common_libraries})	
# add_version(${appname} octree-cpp "" "")
# 
# # ----------------- OPENMP VERSION ----------------- 
# 
add_library(${appname}-openmp  src/cpp/kernels.cpp)
target_link_libraries(${appname}-openmp   ${common_libraries})	
SET_TARGET_PROPERTIES(${appname}-openmp PROPERTIES COMPILE_FLAGS "-fopenmp")
add_version(${appname} openmp "-fopenmp" "-fopenmp")
# 
# # ----------------- OPENMP VERSION OCTREE ----------------- 
# 
# add_library(${appname}-octree-openmp  src/octree/kernels.cpp src/octree/data-structure/octree.hpp)
# target_link_libraries(${appname}-octree-openmp   ${common_libraries})	
# SET_TARGET_PROPERTIES(${appname}-octree-openmp PROPERTIES COMPILE_FLAGS "-fopenmp")
# add_version(${appname} octree-openmp "-fopenmp" "-fopenmp")

# #  ----------------- CUDA VERSION ----------------- 
# if (CUDA_FOUND) 
# 
#   set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++11)
#   set(kfusion_cuda_srcs
#     src/cuda/kernels.cu
#     thirdparty/kfusion.h
#     src/cpp/marching_cubes.cpp
#     )
# 
#   include_directories(${CUDA_INCLUDE_DIRS})
#   cuda_add_library(${appname}-cuda    STATIC  ${kfusion_cuda_srcs} OPTIONS -O3 )
#   target_link_libraries(${appname}-cuda   ${common_libraries})	
#   add_version(${appname} cuda "-DCUDA" "")
# 
# endif(CUDA_FOUND)

#  ----------------- thirdparty ----------------- 

add_subdirectory( thirdparty    )
