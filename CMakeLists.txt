cmake_minimum_required(VERSION 3.21)
project(ct-iprd LANGUAGES C CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)

include(${CMAKE_SOURCE_DIR}/devtools/cmake/Conan.cmake)
run_conan()
include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CT_IPRD_BUILD_TESTS ON CACHE BOOL "Whether to build the c++ unit tests")
else()
    set(CT_IPRD_BUILD_TESTS OFF CACHE BOOL "Whether to build the c++ unit tests")
endif()

set(CT_IPRD_CUDA OFF CACHE BOOL "Whether to add cuda")
if(CT_IPRD_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_EXTENSIONS OFF)
    find_package(CUDAToolkit REQUIRED)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

if(MSVC)
    add_compile_options(/W3 /EHsc /bigobj /permissive-)
endif()

find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(pybind11 REQUIRED)
find_package(tsl-robin-map REQUIRED)

# CUDA::cudart CUDA::cuda_driver  link

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)
target_link_libraries(${PROJECT_NAME} INTERFACE fmt::fmt spdlog::spdlog_header_only tsl::robin_map Threads::Threads ${CMAKE_THREAD_LIBS_INIT})
add_library(ct-iprd::ct-iprd ALIAS ${PROJECT_NAME})

if(CT_IPRD_BUILD_TESTS)
    find_package(Catch2 REQUIRED)
    add_subdirectory(tests)

    add_subdirectory(benchmarks)
endif()

add_subdirectory(src)

add_subdirectory(examples)
