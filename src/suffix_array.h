#ifndef SUFFIX_ARRAY_H
#define SUFFIX_ARRAY_H

#include <string>
#include <vector>

// Function to build the suffix array and LCP array
void build_SA_LCP(const std::string& text, std::vector<int>& suffix_array, std::vector<int>& plcp, std::vector<int>& lcp);

#endif