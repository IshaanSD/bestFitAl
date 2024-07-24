#include <iostream>
#include <vector>
#include <string>
#include <libsais.h>  // Include the libsais header

// Function to build the suffix array using libsais
void build_suffix_array(const std::string& text, std::vector<int>& suffix_array) {
    // Using libsais to generate the suffix array
    libsais(reinterpret_cast<const uint8_t*>(text.c_str()), suffix_array.data(), text.length(), 0, nullptr);
}

// Function to build the LCP array using libsais
void build_plcp_array_with_libsais(const std::string& text, const std::vector<int>& suffix_array, std::vector<int>& plcp_array) {    
    // Using libsais to generate the LCP array
    libsais_plcp(reinterpret_cast<const uint8_t*>(text.c_str()), suffix_array.data(), plcp_array.data(), text.length());
}

// Function to build the LCP array using libsais
void build_lcp_array_with_libsais(const std::string& text, const std::vector<int>& suffix_array, std::vector<int>& plcp_array, std::vector<int>& lcp_array) {    
    // Using libsais to generate the LCP array
    libsais_lcp(plcp_array.data(), suffix_array.data(), lcp_array.data(), text.length());
}

// Helper function to print arrays
template<typename T>
void print_array(const std::vector<T>& array, const std::string& label) {
    std::cout << label << ": ";
    int i = 0;
    for (const auto& val : array) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void build_SA_LCP(const std::string& text, std::vector<int>& suffix_array, std::vector<int>& plcp, std::vector<int>& lcp) {
    suffix_array.resize(text.size());
    plcp.resize(text.size());
    lcp.resize(text.size());
    build_suffix_array(text,suffix_array);
    build_plcp_array_with_libsais(text,suffix_array,plcp);
    build_lcp_array_with_libsais(text,suffix_array,plcp,lcp);
}

// int main() {
//     // Example usage
//     std::string text = "banana";
    
//     std::cout << "Text: " << text << std::endl;
    
//     // Build empty data structures
//     size_t n = text.length();
//     std::vector<int> suffix_array(n);
//     std::vector<int> plcp_array(n);
//     std::vector<int> lcp_array(n);

//     // Build the suffix array
//     build_SA_LCP(text,suffix_array,plcp_array,lcp_array);
//     print_array(suffix_array, "Suffix Array");


//     print_array(lcp_array, "LCP Array (with libsais)");

//     return 0;
// }