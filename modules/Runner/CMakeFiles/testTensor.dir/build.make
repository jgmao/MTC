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
include modules/Runner/CMakeFiles/testTensor.dir/depend.make

# Include the progress variables for this target.
include modules/Runner/CMakeFiles/testTensor.dir/progress.make

# Include the compile flags for this target's objects.
include modules/Runner/CMakeFiles/testTensor.dir/flags.make

modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o: modules/Runner/CMakeFiles/testTensor.dir/flags.make
modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o: modules/Runner/src/testTensor.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/guoxin/Projects/MTC/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o"
	cd /home/guoxin/Projects/MTC/modules/Runner && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/testTensor.dir/src/testTensor.cpp.o -c /home/guoxin/Projects/MTC/modules/Runner/src/testTensor.cpp

modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testTensor.dir/src/testTensor.cpp.i"
	cd /home/guoxin/Projects/MTC/modules/Runner && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/guoxin/Projects/MTC/modules/Runner/src/testTensor.cpp > CMakeFiles/testTensor.dir/src/testTensor.cpp.i

modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testTensor.dir/src/testTensor.cpp.s"
	cd /home/guoxin/Projects/MTC/modules/Runner && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/guoxin/Projects/MTC/modules/Runner/src/testTensor.cpp -o CMakeFiles/testTensor.dir/src/testTensor.cpp.s

modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.requires:
.PHONY : modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.requires

modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.provides: modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.requires
	$(MAKE) -f modules/Runner/CMakeFiles/testTensor.dir/build.make modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.provides.build
.PHONY : modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.provides

modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.provides.build: modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o

# Object files for target testTensor
testTensor_OBJECTS = \
"CMakeFiles/testTensor.dir/src/testTensor.cpp.o"

# External object files for target testTensor
testTensor_EXTERNAL_OBJECTS =

Debug/bin/testTensor: modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o
Debug/bin/testTensor: /usr/local/lib/libopencv_gpuimgproc.so
Debug/bin/testTensor: /usr/local/lib/libopencv_gpuarithm.so
Debug/bin/testTensor: /usr/local/lib/libopencv_gpufilters.so
Debug/bin/testTensor: /usr/local/lib/libopencv_imgproc.so
Debug/bin/testTensor: /usr/local/lib/libopencv_highgui.so
Debug/bin/testTensor: /usr/local/lib/libopencv_core.so
Debug/bin/testTensor: modules/TensorLib/libCube.a
Debug/bin/testTensor: modules/TensorLib/libSize3.a
Debug/bin/testTensor: modules/Runner/libTester.a
Debug/bin/testTensor: modules/TensorLib/libTensorLite.a
Debug/bin/testTensor: modules/Metric/libSteerable.a
Debug/bin/testTensor: modules/Metric/libLRI.a
Debug/bin/testTensor: modules/Metric/libalgorithms.a
Debug/bin/testTensor: Debug/lib/libutil.so
Debug/bin/testTensor: modules/Metric/libLRI.a
Debug/bin/testTensor: modules/TensorLib/libCube.a
Debug/bin/testTensor: modules/TensorLib/libSize3.a
Debug/bin/testTensor: modules/TensorLib/libTensorLite.a
Debug/bin/testTensor: Debug/lib/libutil.so
Debug/bin/testTensor: /usr/local/lib/libopencv_gpuimgproc.so
Debug/bin/testTensor: /usr/local/lib/libopencv_gpuarithm.so
Debug/bin/testTensor: /usr/local/lib/libopencv_gpufilters.so
Debug/bin/testTensor: /usr/local/lib/libopencv_imgproc.so
Debug/bin/testTensor: /usr/local/lib/libopencv_highgui.so
Debug/bin/testTensor: /usr/local/lib/libopencv_core.so
Debug/bin/testTensor: modules/Runner/CMakeFiles/testTensor.dir/build.make
Debug/bin/testTensor: modules/Runner/CMakeFiles/testTensor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../Debug/bin/testTensor"
	cd /home/guoxin/Projects/MTC/modules/Runner && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testTensor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
modules/Runner/CMakeFiles/testTensor.dir/build: Debug/bin/testTensor
.PHONY : modules/Runner/CMakeFiles/testTensor.dir/build

modules/Runner/CMakeFiles/testTensor.dir/requires: modules/Runner/CMakeFiles/testTensor.dir/src/testTensor.cpp.o.requires
.PHONY : modules/Runner/CMakeFiles/testTensor.dir/requires

modules/Runner/CMakeFiles/testTensor.dir/clean:
	cd /home/guoxin/Projects/MTC/modules/Runner && $(CMAKE_COMMAND) -P CMakeFiles/testTensor.dir/cmake_clean.cmake
.PHONY : modules/Runner/CMakeFiles/testTensor.dir/clean

modules/Runner/CMakeFiles/testTensor.dir/depend:
	cd /home/guoxin/Projects/MTC && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/Runner /home/guoxin/Projects/MTC /home/guoxin/Projects/MTC/modules/Runner /home/guoxin/Projects/MTC/modules/Runner/CMakeFiles/testTensor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/Runner/CMakeFiles/testTensor.dir/depend

