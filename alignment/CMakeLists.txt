# Minimum required version of CMake
cmake_minimum_required(VERSION 3.21)

# Project name and language
project(Alignment CXX)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add path to ROOT's CMake modules, replace this with your actual ROOT installation path
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/bin/root")

# Find ROOT and specify the components you need
find_package(ROOT REQUIRED COMPONENTS Core Tree MathCore Hist RIO Graf Gpad)

if (${APPLE})
    find_package( OpenMP )

# Add ROOT includes
include_directories(${ROOT_INCLUDE_DIRS})

# Add ROOT's definitions
add_definitions(${ROOT_DEFINITIONS})

# Define your executable and its sources
add_executable(Alignment alignmentC.cpp)

# Link against ROOT libraries
target_link_libraries(Alignment PUBLIC ${ROOT_LIBRARIES})
