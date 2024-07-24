include("/home/i3gupta/bestFitAl/build/cmake/CPM_0.39.0.cmake")
CPMAddPackage("NAME;libsais;GITHUB_REPOSITORY;IlyaGrebnov/libsais;GIT_TAG;v2.8.1;OPTIONS;LIBSAIS_USE_OPENMP ON;LIBSAIS_BUILD_SHARED_LIB OFF")
set(libsais_FOUND TRUE)