
if(NOT "/home/i3gupta/bestFitAl/build/_deps/libsais-subbuild/libsais-populate-prefix/src/libsais-populate-stamp/libsais-populate-gitinfo.txt" IS_NEWER_THAN "/home/i3gupta/bestFitAl/build/_deps/libsais-subbuild/libsais-populate-prefix/src/libsais-populate-stamp/libsais-populate-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/home/i3gupta/bestFitAl/build/_deps/libsais-subbuild/libsais-populate-prefix/src/libsais-populate-stamp/libsais-populate-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/home/i3gupta/bestFitAl/build/_deps/libsais-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/i3gupta/bestFitAl/build/_deps/libsais-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone --no-checkout --config "advice.detachedHead=false" "https://github.com/IlyaGrebnov/libsais.git" "libsais-src"
    WORKING_DIRECTORY "/home/i3gupta/bestFitAl/build/_deps"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/IlyaGrebnov/libsais.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout v2.8.1 --
  WORKING_DIRECTORY "/home/i3gupta/bestFitAl/build/_deps/libsais-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'v2.8.1'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git"  submodule update --recursive --init 
    WORKING_DIRECTORY "/home/i3gupta/bestFitAl/build/_deps/libsais-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/i3gupta/bestFitAl/build/_deps/libsais-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/home/i3gupta/bestFitAl/build/_deps/libsais-subbuild/libsais-populate-prefix/src/libsais-populate-stamp/libsais-populate-gitinfo.txt"
    "/home/i3gupta/bestFitAl/build/_deps/libsais-subbuild/libsais-populate-prefix/src/libsais-populate-stamp/libsais-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/i3gupta/bestFitAl/build/_deps/libsais-subbuild/libsais-populate-prefix/src/libsais-populate-stamp/libsais-populate-gitclone-lastrun.txt'")
endif()

