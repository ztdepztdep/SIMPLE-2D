# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ztdep/Projects/SimpleCollocatedF

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ztdep/Projects/SimpleCollocatedF/build

# Include any dependencies generated for this target.
include CMakeFiles/simplecollocatedf.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/simplecollocatedf.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simplecollocatedf.dir/flags.make

CMakeFiles/simplecollocatedf.dir/Modules.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Modules.f95.o: ../Modules.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Modules.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Modules.f95 -o CMakeFiles/simplecollocatedf.dir/Modules.f95.o

CMakeFiles/simplecollocatedf.dir/Modules.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Modules.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Modules.f95 > CMakeFiles/simplecollocatedf.dir/Modules.f95.i

CMakeFiles/simplecollocatedf.dir/Modules.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Modules.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Modules.f95 -o CMakeFiles/simplecollocatedf.dir/Modules.f95.s

CMakeFiles/simplecollocatedf.dir/Modules.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Modules.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Modules.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Modules.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Modules.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Modules.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Modules.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Modules.f95.o


CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o: ../Main_Grid.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Main_Grid.f95 -o CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o

CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Main_Grid.f95 > CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.i

CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Main_Grid.f95 -o CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.s

CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o


CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o: ../Main_MIM.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Main_MIM.f95 -o CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o

CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Main_MIM.f95 > CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.i

CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Main_MIM.f95 -o CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.s

CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o


CMakeFiles/simplecollocatedf.dir/Main_P.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Main_P.f95.o: ../Main_P.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Main_P.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Main_P.f95 -o CMakeFiles/simplecollocatedf.dir/Main_P.f95.o

CMakeFiles/simplecollocatedf.dir/Main_P.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Main_P.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Main_P.f95 > CMakeFiles/simplecollocatedf.dir/Main_P.f95.i

CMakeFiles/simplecollocatedf.dir/Main_P.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Main_P.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Main_P.f95 -o CMakeFiles/simplecollocatedf.dir/Main_P.f95.s

CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Main_P.f95.o


CMakeFiles/simplecollocatedf.dir/Main_T.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Main_T.f95.o: ../Main_T.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Main_T.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Main_T.f95 -o CMakeFiles/simplecollocatedf.dir/Main_T.f95.o

CMakeFiles/simplecollocatedf.dir/Main_T.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Main_T.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Main_T.f95 > CMakeFiles/simplecollocatedf.dir/Main_T.f95.i

CMakeFiles/simplecollocatedf.dir/Main_T.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Main_T.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Main_T.f95 -o CMakeFiles/simplecollocatedf.dir/Main_T.f95.s

CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Main_T.f95.o


CMakeFiles/simplecollocatedf.dir/Main_X.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Main_X.f95.o: ../Main_X.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Main_X.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Main_X.f95 -o CMakeFiles/simplecollocatedf.dir/Main_X.f95.o

CMakeFiles/simplecollocatedf.dir/Main_X.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Main_X.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Main_X.f95 > CMakeFiles/simplecollocatedf.dir/Main_X.f95.i

CMakeFiles/simplecollocatedf.dir/Main_X.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Main_X.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Main_X.f95 -o CMakeFiles/simplecollocatedf.dir/Main_X.f95.s

CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Main_X.f95.o


CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o: ../Main_Y.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Main_Y.f95 -o CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o

CMakeFiles/simplecollocatedf.dir/Main_Y.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Main_Y.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Main_Y.f95 > CMakeFiles/simplecollocatedf.dir/Main_Y.f95.i

CMakeFiles/simplecollocatedf.dir/Main_Y.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Main_Y.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Main_Y.f95 -o CMakeFiles/simplecollocatedf.dir/Main_Y.f95.s

CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o


CMakeFiles/simplecollocatedf.dir/Solver.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Solver.f95.o: ../Solver.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Solver.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Solver.f95 -o CMakeFiles/simplecollocatedf.dir/Solver.f95.o

CMakeFiles/simplecollocatedf.dir/Solver.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Solver.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Solver.f95 > CMakeFiles/simplecollocatedf.dir/Solver.f95.i

CMakeFiles/simplecollocatedf.dir/Solver.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Solver.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Solver.f95 -o CMakeFiles/simplecollocatedf.dir/Solver.f95.s

CMakeFiles/simplecollocatedf.dir/Solver.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Solver.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Solver.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Solver.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Solver.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Solver.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Solver.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Solver.f95.o


CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o: ../LIDDRIVEN.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building Fortran object CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/LIDDRIVEN.f95 -o CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o

CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/LIDDRIVEN.f95 > CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.i

CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/LIDDRIVEN.f95 -o CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.s

CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.requires

CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.provides: CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.provides

CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o


CMakeFiles/simplecollocatedf.dir/Main.f95.o: CMakeFiles/simplecollocatedf.dir/flags.make
CMakeFiles/simplecollocatedf.dir/Main.f95.o: ../Main.f95
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building Fortran object CMakeFiles/simplecollocatedf.dir/Main.f95.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/ztdep/Projects/SimpleCollocatedF/Main.f95 -o CMakeFiles/simplecollocatedf.dir/Main.f95.o

CMakeFiles/simplecollocatedf.dir/Main.f95.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/simplecollocatedf.dir/Main.f95.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/ztdep/Projects/SimpleCollocatedF/Main.f95 > CMakeFiles/simplecollocatedf.dir/Main.f95.i

CMakeFiles/simplecollocatedf.dir/Main.f95.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/simplecollocatedf.dir/Main.f95.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/ztdep/Projects/SimpleCollocatedF/Main.f95 -o CMakeFiles/simplecollocatedf.dir/Main.f95.s

CMakeFiles/simplecollocatedf.dir/Main.f95.o.requires:

.PHONY : CMakeFiles/simplecollocatedf.dir/Main.f95.o.requires

CMakeFiles/simplecollocatedf.dir/Main.f95.o.provides: CMakeFiles/simplecollocatedf.dir/Main.f95.o.requires
	$(MAKE) -f CMakeFiles/simplecollocatedf.dir/build.make CMakeFiles/simplecollocatedf.dir/Main.f95.o.provides.build
.PHONY : CMakeFiles/simplecollocatedf.dir/Main.f95.o.provides

CMakeFiles/simplecollocatedf.dir/Main.f95.o.provides.build: CMakeFiles/simplecollocatedf.dir/Main.f95.o


# Object files for target simplecollocatedf
simplecollocatedf_OBJECTS = \
"CMakeFiles/simplecollocatedf.dir/Modules.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Main_P.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Main_T.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Main_X.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Solver.f95.o" \
"CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o" \
"CMakeFiles/simplecollocatedf.dir/Main.f95.o"

# External object files for target simplecollocatedf
simplecollocatedf_EXTERNAL_OBJECTS =

simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Modules.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Main_P.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Main_T.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Main_X.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Solver.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/Main.f95.o
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/build.make
simplecollocatedf: CMakeFiles/simplecollocatedf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking Fortran executable simplecollocatedf"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simplecollocatedf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simplecollocatedf.dir/build: simplecollocatedf

.PHONY : CMakeFiles/simplecollocatedf.dir/build

CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Modules.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Main_Grid.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Main_MIM.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Main_P.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Main_T.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Main_X.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Main_Y.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Solver.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/LIDDRIVEN.f95.o.requires
CMakeFiles/simplecollocatedf.dir/requires: CMakeFiles/simplecollocatedf.dir/Main.f95.o.requires

.PHONY : CMakeFiles/simplecollocatedf.dir/requires

CMakeFiles/simplecollocatedf.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simplecollocatedf.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simplecollocatedf.dir/clean

CMakeFiles/simplecollocatedf.dir/depend:
	cd /home/ztdep/Projects/SimpleCollocatedF/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ztdep/Projects/SimpleCollocatedF /home/ztdep/Projects/SimpleCollocatedF /home/ztdep/Projects/SimpleCollocatedF/build /home/ztdep/Projects/SimpleCollocatedF/build /home/ztdep/Projects/SimpleCollocatedF/build/CMakeFiles/simplecollocatedf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/simplecollocatedf.dir/depend

