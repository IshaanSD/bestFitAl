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
CMAKE_SOURCE_DIR = /home/i3gupta/bestFitAl

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/i3gupta/bestFitAl/build

# Include any dependencies generated for this target.
include _deps/libsais-build/CMakeFiles/libsais.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include _deps/libsais-build/CMakeFiles/libsais.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/libsais-build/CMakeFiles/libsais.dir/progress.make

# Include the compile flags for this target's objects.
include _deps/libsais-build/CMakeFiles/libsais.dir/flags.make

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.o: _deps/libsais-build/CMakeFiles/libsais.dir/flags.make
_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.o: _deps/libsais-src/src/libsais.c
_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.o: _deps/libsais-build/CMakeFiles/libsais.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/i3gupta/bestFitAl/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.o"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.o -MF CMakeFiles/libsais.dir/src/libsais.c.o.d -o CMakeFiles/libsais.dir/src/libsais.c.o -c /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais.c

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/libsais.dir/src/libsais.c.i"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais.c > CMakeFiles/libsais.dir/src/libsais.c.i

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/libsais.dir/src/libsais.c.s"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais.c -o CMakeFiles/libsais.dir/src/libsais.c.s

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.o: _deps/libsais-build/CMakeFiles/libsais.dir/flags.make
_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.o: _deps/libsais-src/src/libsais16.c
_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.o: _deps/libsais-build/CMakeFiles/libsais.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/i3gupta/bestFitAl/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.o"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.o -MF CMakeFiles/libsais.dir/src/libsais16.c.o.d -o CMakeFiles/libsais.dir/src/libsais16.c.o -c /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais16.c

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/libsais.dir/src/libsais16.c.i"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais16.c > CMakeFiles/libsais.dir/src/libsais16.c.i

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/libsais.dir/src/libsais16.c.s"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais16.c -o CMakeFiles/libsais.dir/src/libsais16.c.s

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.o: _deps/libsais-build/CMakeFiles/libsais.dir/flags.make
_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.o: _deps/libsais-src/src/libsais64.c
_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.o: _deps/libsais-build/CMakeFiles/libsais.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/i3gupta/bestFitAl/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.o"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.o -MF CMakeFiles/libsais.dir/src/libsais64.c.o.d -o CMakeFiles/libsais.dir/src/libsais64.c.o -c /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais64.c

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/libsais.dir/src/libsais64.c.i"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais64.c > CMakeFiles/libsais.dir/src/libsais64.c.i

_deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/libsais.dir/src/libsais64.c.s"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && /home/i3gupta/anaconda3/bin/x86_64-conda-linux-gnu-cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/i3gupta/bestFitAl/build/_deps/libsais-src/src/libsais64.c -o CMakeFiles/libsais.dir/src/libsais64.c.s

# Object files for target libsais
libsais_OBJECTS = \
"CMakeFiles/libsais.dir/src/libsais.c.o" \
"CMakeFiles/libsais.dir/src/libsais16.c.o" \
"CMakeFiles/libsais.dir/src/libsais64.c.o"

# External object files for target libsais
libsais_EXTERNAL_OBJECTS =

_deps/libsais-build/liblibsais.a: _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais.c.o
_deps/libsais-build/liblibsais.a: _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais16.c.o
_deps/libsais-build/liblibsais.a: _deps/libsais-build/CMakeFiles/libsais.dir/src/libsais64.c.o
_deps/libsais-build/liblibsais.a: _deps/libsais-build/CMakeFiles/libsais.dir/build.make
_deps/libsais-build/liblibsais.a: _deps/libsais-build/CMakeFiles/libsais.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/i3gupta/bestFitAl/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C static library liblibsais.a"
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && $(CMAKE_COMMAND) -P CMakeFiles/libsais.dir/cmake_clean_target.cmake
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/libsais.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
_deps/libsais-build/CMakeFiles/libsais.dir/build: _deps/libsais-build/liblibsais.a
.PHONY : _deps/libsais-build/CMakeFiles/libsais.dir/build

_deps/libsais-build/CMakeFiles/libsais.dir/clean:
	cd /home/i3gupta/bestFitAl/build/_deps/libsais-build && $(CMAKE_COMMAND) -P CMakeFiles/libsais.dir/cmake_clean.cmake
.PHONY : _deps/libsais-build/CMakeFiles/libsais.dir/clean

_deps/libsais-build/CMakeFiles/libsais.dir/depend:
	cd /home/i3gupta/bestFitAl/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/i3gupta/bestFitAl /home/i3gupta/bestFitAl/build/_deps/libsais-src /home/i3gupta/bestFitAl/build /home/i3gupta/bestFitAl/build/_deps/libsais-build /home/i3gupta/bestFitAl/build/_deps/libsais-build/CMakeFiles/libsais.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : _deps/libsais-build/CMakeFiles/libsais.dir/depend

