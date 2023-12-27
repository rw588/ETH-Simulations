# Minimum required version of CMake
#cmake_minimum_required(VERSION 3.21)

# Specify the compilers (replace with your LLVM installation path)
#set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
#set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

# Project name and language
#project(Alignment CXX)

# Set C++ standard to C++17
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add path to ROOT's CMake modules, replace this with your actual ROOT installation path
#list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/bin/root")

# Find ROOT and specify the components you need
#find_package(ROOT REQUIRED COMPONENTS Core Tree MathCore Hist RIO Graf Gpad)

# Add ROOT includes
#include_directories(${ROOT_INCLUDE_DIRS})

# Add ROOT's definitions
#add_definitions(${ROOT_DEFINITIONS})

# Define your executable and its sources
#add_executable(Alignment alignmentC.cpp)

# Link against ROOT libraries
#target_link_libraries(Alignment PUBLIC ${ROOT_LIBRARIES})

# Check if compiling on Apple system and find OpenMP
#if (${APPLE})
#    find_package(OpenMP)
#    if(OpenMP_CXX_FOUND)
#        target_link_libraries(Alignment PUBLIC OpenMP::OpenMP_CXX)
#    endif()
#endif()


######
cmake_minimum_required(VERSION 3.21)
project(Alignment CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Set C and C++ compilers
set(CMAKE_C_COMPILER "/opt/homebrew/bin/gcc-13")
set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-13")

# Other settings...

# Append ROOT to the CMake prefix path
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/root/6.30.00")

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS Core Tree MathCore Hist RIO Graf Gpad)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Alignment PUBLIC OpenMP::OpenMP_CXX)
endif()

# Include directories from ROOT
include_directories(${ROOT_INCLUDE_DIRS})

# Add executable
add_executable(Alignment alignmentC.cpp)

# Link ROOT libraries
target_link_libraries(Alignment ${ROOT_LIBRARIES})
