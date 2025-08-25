cmake_minimum_required(VERSION 3.14)
project(channel CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../include)

set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/LBMDriver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/lattice/D2Q9.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/lattice/D2Q5.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/lattice/D3Q19.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/collision/BGK.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/boundary/BounceAndBack.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/Geometry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/Streaming.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/Output.cpp
)

add_executable(channel ${SRC})
