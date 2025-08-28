cmake_minimum_required(VERSION 3.14)
project(taylorGreen CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../include)

file(GLOB_RECURSE SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../src/*.cpp")

add_executable(taylorGreen ${SRC})

find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARY NAMES fftw3)

if (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
    message(STATUS "FFTW3 found: ${FFTW3_LIBRARY}")
    target_include_directories(taylorGreen PRIVATE ${FFTW3_INCLUDE_DIR})
    target_link_libraries(taylorGreen PRIVATE ${FFTW3_LIBRARY})
else()
    message(FATAL_ERROR "FFTW3 not found")
endif()

find_package(ZLIB REQUIRED)
target_link_libraries(taylorGreen PRIVATE ZLIB::ZLIB)