# Set the minimum CMake version required for the project
cmake_minimum_required(VERSION 3.10)

# Compiler options
if(compiler STREQUAL intel)
    # Use Intel compiler if specified
    set(CMAKE_C_COMPILER "icx")
    set(CMAKE_CXX_COMPILER "icpx")
    set(base_compile_options "-Wall -Wextra -O3")
elseif(DEFINED compiler)
    # Error if an invalid compiler is specified
    message(FATAL ERROR "invalid compiler")
endif()

# Project name
project(HelloWorld)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)

# Source files
file(GLOB SOURCES
    "src/main.cpp"
    "src/Greeter.cpp"
)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add executable
add_executable(HelloWorld ${SOURCES})