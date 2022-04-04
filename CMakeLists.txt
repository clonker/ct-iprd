cmake_minimum_required(VERSION 3.9)
project(ct-iprd LANGUAGES C CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)

set(CT_IPRD_BUILD_TESTS OFF CACHE BOOL "Whether to build the c++ unit tests")

include(${CMAKE_SOURCE_DIR}/devtools/cmake/Conan.cmake)
run_conan()
include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

if(MSVC)
    add_compile_options(/W3 /EHsc /bigobj /permissive-)
endif()

find_package(pybind11 REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(Threads REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)
target_link_libraries(${PROJECT_NAME} INTERFACE fmt::fmt spdlog::spdlog_header_only Threads::Threads ${CMAKE_THREAD_LIBS_INIT})
add_library(ct-iprd::ct-iprd ALIAS ${PROJECT_NAME})

if(CT_IPRD_BUILD_TESTS)
    find_package(Catch2 REQUIRED)
    add_subdirectory(tests)
endif()

add_subdirectory(examples)
