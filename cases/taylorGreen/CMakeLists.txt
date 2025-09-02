cmake_minimum_required(VERSION 3.14)

set(PROBLEM_NAME lbm)

project(${PROBLEM_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

file(GLOB_RECURSE SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../src/*.cpp")

add_executable(${PROBLEM_NAME} ${SRC})

target_include_directories(${PROBLEM_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
)

# --- FFTW3 ---
find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARY NAMES fftw3)
if (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
    target_include_directories(${PROBLEM_NAME} PRIVATE ${FFTW3_INCLUDE_DIR})
    target_link_libraries(${PROBLEM_NAME} PRIVATE ${FFTW3_LIBRARY})
    message(STATUS "FFTW3 found: ${FFTW3_LIBRARY}")
else()
    message(FATAL_ERROR "FFTW3 not found")
endif()

# --- ZLIB ---
find_package(ZLIB REQUIRED)
target_link_libraries(${PROBLEM_NAME} PRIVATE ZLIB::ZLIB)

# --- OpenMP ---
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROBLEM_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()