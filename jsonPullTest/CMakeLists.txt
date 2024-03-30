cmake_minimum_required(VERSION 3.0)
project(jsonTest)

# Find the nlohmann json library
find_package(nlohmann_json REQUIRED)

# Include directories
include_directories(${nlohmann_json_INCLUDE_DIRS})

# Add the main source file
add_executable(jsonTest jsonTest.cpp)

# Link the nlohmann json library
target_link_libraries(jsonTest ${nlohmann_json_LIBRARIES})
