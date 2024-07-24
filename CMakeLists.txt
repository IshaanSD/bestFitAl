cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)



set(CMAKE_CXX_FLAGS "-lstdc++fs -ggdb3" )

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc static-libstdc++ -Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -ltbb -fsigned-char -fsanitize=address")

project(BestFitAlProject)

# Add libsais
include(CPM.cmake)
CPMAddPackage(
  NAME libsais
  GITHUB_REPOSITORY IlyaGrebnov/libsais
  GIT_TAG v2.8.1
  OPTIONS
    "LIBSAIS_USE_OPENMP ON"
    "LIBSAIS_BUILD_SHARED_LIB OFF"
)

add_subdirectory(src)