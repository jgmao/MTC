# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/guoxin/Projects/MTC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/guoxin/Projects/MTC

# Include any dependencies generated for this target.
include modules/TensorLib/CMakeFiles/Tensor.dir/depend.make

# Include the progress variables for this target.
include modules/TensorLib/CMakeFiles/Tensor.dir/progress.make

# Include the compile flags for this target's objects.
include modules/TensorLib/CMakeFiles/Tensor.dir/flags.make

modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o: modules/TensorLib/CMakeFiles/Tensor.dir/flags.make
modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o: modules/TensorLib/src/Tensor.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/guoxin/Projects/MTC/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Tensor.dir/src/Tensor.cpp.o -c /home/guoxin/Projects/MTC/modules/TensorLib/src/Tensor.cpp

modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tensor.dir/src/Tensor.cpp.i"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/guoxin/Projects/MTC/modules/TensorLib/src/Tensor.cpp > CMakeFiles/Tensor.dir/src/Tensor.cpp.i

modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tensor.dir/src/Tensor.cpp.s"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/guoxin/Projects/MTC/modules/TensorLib/src/Tensor.cpp -o CMakeFiles/Tensor.dir/src/Tensor.cpp.s

modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.requires:
.PHONY : modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.requires

modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.provides: modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.requires
	$(MAKE) -f modules/TensorLib/CMakeFiles/Tensor.dir/build.make modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.provides.build
.PHONY : modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.provides

modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.provides.build: modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o

# Object files for target Tensor
Tensor_OBJECTS = \
"CMakeFiles/Tensor.dir/src/Tensor.cpp.o"

# External object files for target Tensor
Tensor_EXTERNAL_OBJECTS =

modules/TensorLib/libTensor.a: modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o
modules/TensorLib/libTensor.a: modules/TensorLib/CMakeFiles/Tensor.dir/build.make
modules/TensorLib/libTensor.a: modules/TensorLib/CMakeFiles/Tensor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libTensor.a"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && $(CMAKE_COMMAND) -P CMakeFiles/Tensor.dir/cmake_clean_target.cmake
	cd /home/guoxin/Projects/MTC/modules/TensorLib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Tensor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
modules/TensorLib/CMakeFiles/Tensor.dir/build: modules/TensorLib/libTensor.a
.PHONY : modules/TensorLib/CMakeFiles/Tensor.dir/build

modules/TensorLib/CMakeFiles/Tensor.dir/requires: modules/TensorLib/CMakeFiles/Tensor.dir/src/Tensor.cpp.o.requires
.PHONY : modules/TensorLib/CMakeFiles/Tensor.dir/requires

modules/TensorLib/CMakeFiles/Tensor.dir/clean:
	cd /home/guoxin/Projects/MTC/modules/TensorLib && $(CMAKE_COMMAND) -P CMakeFiles/Tensor.dir/cmake_clean.cmake
.PHONY : modules/TensorLib/CMakeFiles/Tensor.dir/clean

modules/TensorLib/CMakeFiles/Tensor.dir/depend:
	cd /home/guoxin/Projects/MTC && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/TensorLib /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/TensorLib /home/guoxin/Projects/MTC/modules/TensorLib/CMakeFiles/Tensor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/TensorLib/CMakeFiles/Tensor.dir/depend
