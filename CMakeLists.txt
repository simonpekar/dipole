cmake_minimum_required(VERSION 2.8)
project(P3A)

# Use C++11
set (CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT COMPONENTS MathMore REQUIRED)
include(${ROOT_USE_FILE})

#find_package(GSL REQUIRED)
#include_directories(${GSL_INCLUDE_DIRS})

add_library(build build.cpp)
target_link_libraries(build ${ROOT_LIBRARIES})
target_include_directories(build PUBLIC)

add_library(analyse analyse.cpp)
target_link_libraries(analyse ${ROOT_LIBRARIES})
target_include_directories(analyse PUBLIC)

add_executable(main main.cpp)
target_link_libraries(main ${ROOT_LIBRARIES} build)
target_include_directories(main PUBLIC)