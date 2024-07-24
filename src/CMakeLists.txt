set(BFA_FILES
    best_fit_al.cpp
    suffix_array.cpp
)
# Define best_fit_al target
add_executable(BestFitExecutable ${BFA_FILES})
target_link_libraries(BestFitExecutable libsais)

# Other targets

# sufix_array
add_executable(SuffixArrayExecutable suffix_array.cpp)
target_link_libraries(SuffixArrayExecutable libsais)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})