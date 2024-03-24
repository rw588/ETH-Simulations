# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rob/Documents/Coding/ETH-Simulations/alignment

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rob/Documents/Coding/ETH-Simulations/alignment/build

# Include any dependencies generated for this target.
include CMakeFiles/Alignment.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Alignment.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Alignment.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Alignment.dir/flags.make

CMakeFiles/Alignment.dir/alignmentC.cpp.o: CMakeFiles/Alignment.dir/flags.make
CMakeFiles/Alignment.dir/alignmentC.cpp.o: ../alignmentC.cpp
CMakeFiles/Alignment.dir/alignmentC.cpp.o: CMakeFiles/Alignment.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rob/Documents/Coding/ETH-Simulations/alignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Alignment.dir/alignmentC.cpp.o"
	/home/rob/miniconda3/envs/pyROOT/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Alignment.dir/alignmentC.cpp.o -MF CMakeFiles/Alignment.dir/alignmentC.cpp.o.d -o CMakeFiles/Alignment.dir/alignmentC.cpp.o -c /home/rob/Documents/Coding/ETH-Simulations/alignment/alignmentC.cpp

CMakeFiles/Alignment.dir/alignmentC.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Alignment.dir/alignmentC.cpp.i"
	/home/rob/miniconda3/envs/pyROOT/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rob/Documents/Coding/ETH-Simulations/alignment/alignmentC.cpp > CMakeFiles/Alignment.dir/alignmentC.cpp.i

CMakeFiles/Alignment.dir/alignmentC.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Alignment.dir/alignmentC.cpp.s"
	/home/rob/miniconda3/envs/pyROOT/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rob/Documents/Coding/ETH-Simulations/alignment/alignmentC.cpp -o CMakeFiles/Alignment.dir/alignmentC.cpp.s

# Object files for target Alignment
Alignment_OBJECTS = \
"CMakeFiles/Alignment.dir/alignmentC.cpp.o"

# External object files for target Alignment
Alignment_EXTERNAL_OBJECTS =

Alignment: CMakeFiles/Alignment.dir/alignmentC.cpp.o
Alignment: CMakeFiles/Alignment.dir/build.make
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libCore.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libImt.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libRIO.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libNet.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libHist.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libGraf.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libGraf3d.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libGpad.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libROOTDataFrame.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libTree.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libTreePlayer.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libRint.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libPostscript.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libMatrix.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libPhysics.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libMathCore.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libThread.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libMultiProc.so
Alignment: /home/rob/miniconda3/envs/pyROOT/lib/libROOTVecOps.so
Alignment: CMakeFiles/Alignment.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rob/Documents/Coding/ETH-Simulations/alignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Alignment"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Alignment.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Alignment.dir/build: Alignment
.PHONY : CMakeFiles/Alignment.dir/build

CMakeFiles/Alignment.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Alignment.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Alignment.dir/clean

CMakeFiles/Alignment.dir/depend:
	cd /home/rob/Documents/Coding/ETH-Simulations/alignment/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rob/Documents/Coding/ETH-Simulations/alignment /home/rob/Documents/Coding/ETH-Simulations/alignment /home/rob/Documents/Coding/ETH-Simulations/alignment/build /home/rob/Documents/Coding/ETH-Simulations/alignment/build /home/rob/Documents/Coding/ETH-Simulations/alignment/build/CMakeFiles/Alignment.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Alignment.dir/depend

