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
include modules/TensorLib/CMakeFiles/Size3.dir/depend.make

# Include the progress variables for this target.
include modules/TensorLib/CMakeFiles/Size3.dir/progress.make

# Include the compile flags for this target's objects.
include modules/TensorLib/CMakeFiles/Size3.dir/flags.make

modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o: modules/TensorLib/CMakeFiles/Size3.dir/flags.make
modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o: modules/TensorLib/src/Size3.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/guoxin/Projects/MTC/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Size3.dir/src/Size3.cpp.o -c /home/guoxin/Projects/MTC/modules/TensorLib/src/Size3.cpp

modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Size3.dir/src/Size3.cpp.i"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/guoxin/Projects/MTC/modules/TensorLib/src/Size3.cpp > CMakeFiles/Size3.dir/src/Size3.cpp.i

modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Size3.dir/src/Size3.cpp.s"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/guoxin/Projects/MTC/modules/TensorLib/src/Size3.cpp -o CMakeFiles/Size3.dir/src/Size3.cpp.s

modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.requires:
.PHONY : modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.requires

modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.provides: modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.requires
	$(MAKE) -f modules/TensorLib/CMakeFiles/Size3.dir/build.make modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.provides.build
.PHONY : modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.provides

modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.provides.build: modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o

# Object files for target Size3
Size3_OBJECTS = \
"CMakeFiles/Size3.dir/src/Size3.cpp.o"

# External object files for target Size3
Size3_EXTERNAL_OBJECTS =

modules/TensorLib/libSize3.a: modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o
modules/TensorLib/libSize3.a: modules/TensorLib/CMakeFiles/Size3.dir/build.make
modules/TensorLib/libSize3.a: modules/TensorLib/CMakeFiles/Size3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libSize3.a"
	cd /home/guoxin/Projects/MTC/modules/TensorLib && $(CMAKE_COMMAND) -P CMakeFiles/Size3.dir/cmake_clean_target.cmake
	cd /home/guoxin/Projects/MTC/modules/TensorLib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Size3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
modules/TensorLib/CMakeFiles/Size3.dir/build: modules/TensorLib/libSize3.a
.PHONY : modules/TensorLib/CMakeFiles/Size3.dir/build

modules/TensorLib/CMakeFiles/Size3.dir/requires: modules/TensorLib/CMakeFiles/Size3.dir/src/Size3.cpp.o.requires
.PHONY : modules/TensorLib/CMakeFiles/Size3.dir/requires

modules/TensorLib/CMakeFiles/Size3.dir/clean:
	cd /home/guoxin/Projects/MTC/modules/TensorLib && $(CMAKE_COMMAND) -P CMakeFiles/Size3.dir/cmake_clean.cmake
.PHONY : modules/TensorLib/CMakeFiles/Size3.dir/clean

modules/TensorLib/CMakeFiles/Size3.dir/depend:
	cd /home/guoxin/Projects/MTC && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/TensorLib /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/TensorLib /home/guoxin/Projects/MTC/modules/TensorLib/CMakeFiles/Size3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/TensorLib/CMakeFiles/Size3.dir/depend

