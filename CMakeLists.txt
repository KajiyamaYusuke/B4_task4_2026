cmake_minimum_required(VERSION 3.10)

if(compiler STREQUAL intel)
    set(CMAKE_C_COMPILER "icx")
    set(CMAKE_CXX_COMPILER "icpx")
    set(base_compile_options "-Wall -Wextra -O3")
elseif(DEFINED compiler)
    message(FATAL_ERROR "invalid compiler")
endif()

project(LaplaceFEM2D CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(LaplaceFEM2D
    src/laplace.cpp
)
