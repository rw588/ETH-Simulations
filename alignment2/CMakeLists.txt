# Minimum required version of CMake
cmake_minimum_required(VERSION 3.21)

# Project name and language
project(Alignment2 CXX)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add path to ROOT's CMake modules, replace this with your actual ROOT installation path
list(APPEND CMAKE_PREFIX_PATH "/home/rob/miniconda3/envs/pyROOT/bin/root")

# Find ROOT and specify the components you need
find_package(ROOT REQUIRED COMPONENTS Core Tree MathCore Hist RIO Graf Gpad)

# Find the nlohmann json library
find_package(nlohmann_json REQUIRED)

# Define your executable and its sources
add_executable(Alignment2 alignment2.cpp)

# Include ROOT's and nlohmann_json's include directories
target_include_directories(Alignment2 PUBLIC ${ROOT_INCLUDE_DIRS} ${nlohmann_json_INCLUDE_DIRS})

# Link against ROOT libraries and nlohmann_json
target_link_libraries(Alignment2 PUBLIC ${ROOT_LIBRARIES} ${nlohmann_json_LIBRARIES})
