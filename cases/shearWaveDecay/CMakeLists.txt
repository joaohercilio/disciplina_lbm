cmake_minimum_required(VERSION 3.14)
project(shearWaveDecay CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

get_filename_component(LBM_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)

include_directories(${LBM_PROJECT_ROOT}/include)
include_directories(${LBM_PROJECT_ROOT})

set(COMMON_SRC
    ${LBM_PROJECT_ROOT}/src/main.cpp
    ${LBM_PROJECT_ROOT}/src/LBMDriver.cpp
    ${LBM_PROJECT_ROOT}/src/lattice/D2Q9.cpp
    ${LBM_PROJECT_ROOT}/src/lattice/D2Q5.cpp
    ${LBM_PROJECT_ROOT}/src/collision/BGK.cpp
    ${LBM_PROJECT_ROOT}/src/Geometry.cpp
    ${LBM_PROJECT_ROOT}/src/Streaming.cpp
    ${LBM_PROJECT_ROOT}/src/Output.cpp
)

add_executable(shearWaveDecay ${COMMON_SRC})

target_include_directories(shearWaveDecay PRIVATE ${CMAKE_CURRENT_LIST_DIR})
