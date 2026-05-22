cmake_minimum_required(VERSION 3.10)

project(FEMSolver)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler options
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -O2)
endif()

# Eigen
find_package(Eigen3 CONFIG REQUIRED)

# Include directories
include_directories(include)

# Source files
file(GLOB SOURCES
    "src/main.cpp"
    "src/Mesh.cpp"
    "src/FEMSolver.cpp"
)

# Executable
add_executable(FEMSolver ${SOURCES})

# Link Eigen
target_link_libraries(FEMSolver PRIVATE Eigen3::Eigen)