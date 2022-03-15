cmake_minimum_required(VERSION 3.9)
project(ct-iprd LANGUAGES C CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)

set(CT_IPRD_BUILD_TESTS OFF CACHE BOOL "Whether to build the c++ unit tests")


find_package(OpenMP 4)
if(OpenMP_FOUND)
    add_definitions(-DUSE_OPENMP)
endif()

if(MSVC)
    add_compile_options(/W3 /EHsc /bigobj /permissive- /std:c++17)
endif()

find_package(pybind11 REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)
target_link_libraries(${PROJECT_NAME} INTERFACE fmt::fmt spdlog::spdlog)
add_library(ct-iprd::ct-iprd ALIAS ${PROJECT_NAME})

if(CT_IPRD_BUILD_TESTS)
    add_subdirectory(tests)
endif()

add_subdirectory(examples)
