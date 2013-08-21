# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/cmake-gui -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/guoxin/Projects/MTC/CMakeFiles /home/guoxin/Projects/MTC/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/guoxin/Projects/MTC/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named Cube

# Build rule for target.
Cube: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 Cube
.PHONY : Cube

# fast build rule for target.
Cube/fast:
	$(MAKE) -f modules/TensorLib/CMakeFiles/Cube.dir/build.make modules/TensorLib/CMakeFiles/Cube.dir/build
.PHONY : Cube/fast

#=============================================================================
# Target rules for targets named Size3

# Build rule for target.
Size3: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 Size3
.PHONY : Size3

# fast build rule for target.
Size3/fast:
	$(MAKE) -f modules/TensorLib/CMakeFiles/Size3.dir/build.make modules/TensorLib/CMakeFiles/Size3.dir/build
.PHONY : Size3/fast

#=============================================================================
# Target rules for targets named TensorLite

# Build rule for target.
TensorLite: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 TensorLite
.PHONY : TensorLite

# fast build rule for target.
TensorLite/fast:
	$(MAKE) -f modules/TensorLib/CMakeFiles/TensorLite.dir/build.make modules/TensorLib/CMakeFiles/TensorLite.dir/build
.PHONY : TensorLite/fast

#=============================================================================
# Target rules for targets named util

# Build rule for target.
util: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 util
.PHONY : util

# fast build rule for target.
util/fast:
	$(MAKE) -f modules/Utility/CMakeFiles/util.dir/build.make modules/Utility/CMakeFiles/util.dir/build
.PHONY : util/fast

#=============================================================================
# Target rules for targets named LRI

# Build rule for target.
LRI: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 LRI
.PHONY : LRI

# fast build rule for target.
LRI/fast:
	$(MAKE) -f modules/Metric/CMakeFiles/LRI.dir/build.make modules/Metric/CMakeFiles/LRI.dir/build
.PHONY : LRI/fast

#=============================================================================
# Target rules for targets named Steerable

# Build rule for target.
Steerable: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 Steerable
.PHONY : Steerable

# fast build rule for target.
Steerable/fast:
	$(MAKE) -f modules/Metric/CMakeFiles/Steerable.dir/build.make modules/Metric/CMakeFiles/Steerable.dir/build
.PHONY : Steerable/fast

#=============================================================================
# Target rules for targets named algorithms

# Build rule for target.
algorithms: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 algorithms
.PHONY : algorithms

# fast build rule for target.
algorithms/fast:
	$(MAKE) -f modules/Metric/CMakeFiles/algorithms.dir/build.make modules/Metric/CMakeFiles/algorithms.dir/build
.PHONY : algorithms/fast

#=============================================================================
# Target rules for targets named Tester

# Build rule for target.
Tester: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 Tester
.PHONY : Tester

# fast build rule for target.
Tester/fast:
	$(MAKE) -f modules/Runner/CMakeFiles/Tester.dir/build.make modules/Runner/CMakeFiles/Tester.dir/build
.PHONY : Tester/fast

#=============================================================================
# Target rules for targets named gpussim

# Build rule for target.
gpussim: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 gpussim
.PHONY : gpussim

# fast build rule for target.
gpussim/fast:
	$(MAKE) -f modules/Runner/CMakeFiles/gpussim.dir/build.make modules/Runner/CMakeFiles/gpussim.dir/build
.PHONY : gpussim/fast

#=============================================================================
# Target rules for targets named main

# Build rule for target.
main: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 main
.PHONY : main

# fast build rule for target.
main/fast:
	$(MAKE) -f modules/Runner/CMakeFiles/main.dir/build.make modules/Runner/CMakeFiles/main.dir/build
.PHONY : main/fast

#=============================================================================
# Target rules for targets named testTensor

# Build rule for target.
testTensor: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 testTensor
.PHONY : testTensor

# fast build rule for target.
testTensor/fast:
	$(MAKE) -f modules/Runner/CMakeFiles/testTensor.dir/build.make modules/Runner/CMakeFiles/testTensor.dir/build
.PHONY : testTensor/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... Cube"
	@echo "... Size3"
	@echo "... TensorLite"
	@echo "... util"
	@echo "... LRI"
	@echo "... Steerable"
	@echo "... algorithms"
	@echo "... Tester"
	@echo "... gpussim"
	@echo "... main"
	@echo "... testTensor"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

