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
include modules/Runner/CMakeFiles/gpussim.dir/depend.make

# Include the progress variables for this target.
include modules/Runner/CMakeFiles/gpussim.dir/progress.make

# Include the compile flags for this target's objects.
include modules/Runner/CMakeFiles/gpussim.dir/flags.make

modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o: modules/Runner/CMakeFiles/gpussim.dir/flags.make
modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o: modules/Runner/src/gpussim.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/guoxin/Projects/MTC/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o"
	cd /home/guoxin/Projects/MTC/modules/Runner && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/gpussim.dir/src/gpussim.cpp.o -c /home/guoxin/Projects/MTC/modules/Runner/src/gpussim.cpp

modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gpussim.dir/src/gpussim.cpp.i"
	cd /home/guoxin/Projects/MTC/modules/Runner && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/guoxin/Projects/MTC/modules/Runner/src/gpussim.cpp > CMakeFiles/gpussim.dir/src/gpussim.cpp.i

modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gpussim.dir/src/gpussim.cpp.s"
	cd /home/guoxin/Projects/MTC/modules/Runner && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/guoxin/Projects/MTC/modules/Runner/src/gpussim.cpp -o CMakeFiles/gpussim.dir/src/gpussim.cpp.s

modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.requires:
.PHONY : modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.requires

modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.provides: modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.requires
	$(MAKE) -f modules/Runner/CMakeFiles/gpussim.dir/build.make modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.provides.build
.PHONY : modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.provides

modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.provides.build: modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o

# Object files for target gpussim
gpussim_OBJECTS = \
"CMakeFiles/gpussim.dir/src/gpussim.cpp.o"

# External object files for target gpussim
gpussim_EXTERNAL_OBJECTS =

Debug/bin/gpussim: modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o
Debug/bin/gpussim: /usr/local/lib/libopencv_gpuarithm.so
Debug/bin/gpussim: /usr/local/lib/libopencv_gpufilters.so
Debug/bin/gpussim: /usr/local/lib/libopencv_imgproc.so
Debug/bin/gpussim: /usr/local/lib/libopencv_highgui.so
Debug/bin/gpussim: /usr/local/lib/libopencv_core.so
Debug/bin/gpussim: modules/Runner/CMakeFiles/gpussim.dir/build.make
Debug/bin/gpussim: modules/Runner/CMakeFiles/gpussim.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../Debug/bin/gpussim"
	cd /home/guoxin/Projects/MTC/modules/Runner && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gpussim.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
modules/Runner/CMakeFiles/gpussim.dir/build: Debug/bin/gpussim
.PHONY : modules/Runner/CMakeFiles/gpussim.dir/build

modules/Runner/CMakeFiles/gpussim.dir/requires: modules/Runner/CMakeFiles/gpussim.dir/src/gpussim.cpp.o.requires
.PHONY : modules/Runner/CMakeFiles/gpussim.dir/requires

modules/Runner/CMakeFiles/gpussim.dir/clean:
	cd /home/guoxin/Projects/MTC/modules/Runner && $(CMAKE_COMMAND) -P CMakeFiles/gpussim.dir/cmake_clean.cmake
.PHONY : modules/Runner/CMakeFiles/gpussim.dir/clean

modules/Runner/CMakeFiles/gpussim.dir/depend:
	cd /home/guoxin/Projects/MTC && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/Runner /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/Runner /home/guoxin/Projects/MTC/modules/Runner/CMakeFiles/gpussim.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/Runner/CMakeFiles/gpussim.dir/depend

