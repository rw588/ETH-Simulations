cmake_minimum_required(VERSION 3.21)
project(Alignment)

set(CMAKE_CXX_STANDARD 17)

#linux
set(ROOT_LIBRARIES "/home/rob/miniconda3/envs/pyROOT") # Replace with the actual paths

# Use libc++ on macOS
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

# Find ROOT and specify the necessary components
#for mac
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#for linux
list(APPEND CMAKE_PREFIX_PATH /home/rob/miniconda3/envs/pyROOT)

find_package(ROOT REQUIRED COMPONENTS Core Tree MathCore Hist RIO Graf Gpad) # Add more components if needed

# Include ROOT directories and definitions
include_directories(${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})
add_definitions(${ROOT_DEFINITIONS})

# Define the executable and its source files
add_executable(Alignment alignmentC.cpp)

# Include additional directories (if they exist in your project)
target_include_directories(Alignment PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link ROOT libraries, the C++ library, and any other necessary libraries
#mac
#target_link_libraries(Alignment PUBLIC ${ROOT_LIBRARIES} c++)
#linux
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
target_link_libraries(Alignment PUBLIC ${ROOT_LIBRARIES})

# Check for OpenMP and add flags if found
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
