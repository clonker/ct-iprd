cmake_minimum_required(VERSION 3.16)
project(ct-iprd LANGUAGES C CXX VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 20)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CT_IPRD_BUILD_TESTS ON CACHE BOOL "Whether to build the c++ unit tests")
else()
    set(CT_IPRD_BUILD_TESTS OFF CACHE BOOL "Whether to build the c++ unit tests")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
include(${CMAKE_SOURCE_DIR}/devtools/cmake/Conan.cmake)
run_conan()
include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

# Constrain FindPython to find the Python version used by scikit-build
if (SKBUILD)
    message(STATUS "Python_VERSION ${PYTHON_VERSION_STRING}")
    message(STATUS "Python_EXECUTABLE ${PYTHON_EXECUTABLE}")
    message(STATUS "Python_INCLUDE_DIR ${PYTHON_INCLUDE_DIR}")
    message(STATUS "Python_LIBRARIES ${PYTHON_LIBRARY}")
    set(Python_VERSION "${PYTHON_VERSION_STRING}")
    set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}")
    set(Python_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}")
    set(Python_LIBRARIES "${PYTHON_LIBRARY}")
endif()
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

if(MSVC)
    add_compile_options(/W3 /EHsc /bigobj /permissive-)
endif()

add_subdirectory(contrib/nanobind)

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
