#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include "suffix_array.h"  // Include the header where build_SA_LCP is declared

using namespace std;

int max3(int a, int b, int c) {
    return max(max(a, b), c);
}
void write2DArrayToCSV(const std::vector<std::vector<int>>& array, const std::string& filename) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < array.size(); ++i) {
            for (int j = 0; j < array[i].size(); ++j) {
                file << array[i][j];
                if (j < array[i].size() - 1) {
                    file << ",";  // Separate elements in a row with commas
                }
            }
            file << "\n";  // End of row
        }
        file.close();
        std::cout << "2D array successfully written to " << filename << std::endl;
    } else {
        std::cerr << "Failed to open file " << filename << " for writing." << std::endl;
    }
}
pair<int, int> find_10mer_hits(const string& text, const int* suffix_array, const int* lcp, size_t len1, size_t len2) {
    int K = 10;
    int B = 100;
    int gap = -1*B/K;
    int match = 20;
    int mismatch = 0;
    size_t n = text.size();
    map<int, vector<int>> hits_map; // Map to store indices of runs of suffixes sharing 10-mers

    // Iterate over the suffix array to find runs of suffixes sharing 10-mers
    {
    int i = 1;
    int run_start = 0;
    while (i < n) {
        if (lcp[i] < K){
            run_start = i;
        }
        if (lcp[i] >= K) { // Check if the LCP value is greater than or equal to 10
            if (run_start == i-1) hits_map[suffix_array[run_start]].push_back(suffix_array[i-1]);
            hits_map[suffix_array[run_start]].push_back(suffix_array[i]);
        }
        i++;
    }
    }
    map<int, vector<int>> valid_hits_map1;
    map<int, vector<int>> valid_hits_map2;

    // for (const auto& pair : hits_map) {
    //     // print
    //     cout << "Run1 starting at index " << pair.first << ": ";
    //     for (int index : pair.second) {
    //     // print
    //         cout << index << " ";
    //     }
    //     cout << endl;
    // }
    for (const auto& pair : hits_map) {
        vector<int>ind1;
        vector<int>ind2;
        for (int index : pair.second) {

            if (index < len1){
                ind1.push_back(index);
            }
            if (index > len1){
                ind2.push_back(index - len1 - 1);
            }
        }
        if (ind1.size()>0 && ind2.size()>0){
            for (int x : ind1) valid_hits_map1[suffix_array[pair.first]].push_back(x);
            for (int x : ind2) valid_hits_map2[suffix_array[pair.first]].push_back(x);
        }
    }
    // for (const auto& pair : valid_hits_map1) {
    //     // print
    //     cout << "Run1 starting at index " << pair.first << ": ";
    //     for (int index : pair.second) {
    //     // print
    //         cout << index << " ";
    //     }
    //     cout << endl;

    //     // print
    //     cout << "Run2 starting at index " << pair.first << ": ";
    //     for (int index : valid_hits_map2.at(pair.first)) {
    //     // print
    //         cout << index << " ";
    //     }
    //     cout << endl;
    // }

    vector<int> diag_kmer_counts((len1 + len2) / B + 1, 0);
    // Approach 1: Diagonal window-based
    for (const auto& pair : valid_hits_map1) {
        for (int i : pair.second) {
            for (int j : valid_hits_map2.at(pair.first)) {
                int diag = (j-i + len1) / B;
                diag_kmer_counts[diag]++;
            }
        }
    }
    int max_diag=0;
    int max_diag_count = diag_kmer_counts[0];

    // std::ofstream diag_file("diag_scores.csv");
    for (int i = 0; i< diag_kmer_counts.size(); i++){
        if (max_diag_count< diag_kmer_counts[i]){
            max_diag=i;
            max_diag_count = diag_kmer_counts[i];
        }

        // if (i>0) diag_file<<',';
        // diag_file<<diag_kmer_counts[i];
    }
    // diag_file<<endl;
    int start2 = max((int)(max_diag*B - len1),0);
    int end2 = min((int)len2,max_diag*B);
    int end1 = (max_diag*B <= len2) ? len1 : len1 - (max_diag*B - len2);
    int start1 = (max_diag*B > len1) ? 0 : (len1 - max_diag*B);
    cout << start1 << "\n";
    cout << end1 << "\n";
    cout << start2 << "\n";
    cout << end2 << "\n";
    // // Approach 2: Bin-based
    // int nb1 = len1 / B + 1;
    // int nb2 = len2 / B + 1;
    // vector<vector<int>> kmer_counts(nb1, vector<int>(nb2, 0));
    // for (const auto& pair : valid_hits_map1) {
    //     for (int i : pair.second) {
    //         for (int j : valid_hits_map2.at(pair.first)) {
    //             int bin_i = i / B;
    //             int bin_j = j / B;
    //             kmer_counts[bin_i][bin_j]++;
    //             cout << kmer_counts[bin_i][bin_j];
    //         }
    //     }
    // }

    // cout<<"kmer"<<kmer_counts[1][459]<<endl;
    // // fitting al
    // vector<vector<int>> score(nb1, vector<int>(nb2, 0));
    // vector<vector<int>> bt(nb1, vector<int>(nb2, 0));

    // for (int i = 1; i < nb1; ++i) {
    //     score[i][0] = score[i - 1][0] + gap;
    // }
    // // Fill the score matrix
    // for (int i = 1; i < nb1; ++i) {
    //     for (int j = 1; j < nb2; ++j) {
    //         int matchMismatch = score[i - 1][j - 1] + (match*kmer_counts[i][j] + mismatch*max( 0,(B - kmer_counts[i][j]) ));
    //         int deleteFromA = score[i - 1][j] + gap;
    //         int insertToA = score[i][j - 1] + gap;
    //         if (score[i][j - 1] == 0) insertToA = max(0,insertToA);
    //         score[i][j] = kmer_counts[i][j]==0 ? max(deleteFromA, insertToA) : max3(matchMismatch, deleteFromA, insertToA);

    //         // Track the direction from which the score was derived
    //         if (score[i][j] == matchMismatch) {
    //             bt[i][j] = 0; // Diagonal
    //         } else if (score[i][j] == insertToA) {
    //             bt[i][j] = 1; // left
    //         } else {
    //             bt[i][j] = 2; // up
    //         }
    //     }
    // }
    // // backtrack
    // int maxJ = 0;
    // int maxScore = score[nb1-1][0];    
    // for (int j = 1; j < nb2; ++j) {
    //     if (score[nb1-1][j] > maxScore) {
    //         maxScore = score[nb1-1][j];
    //         maxJ = j;
    //     }
    // }

    // int i = nb1-1, j = maxJ;
    // while (i > 0 && j > 0) {
    //     if (bt[i][j] == 0) {
    //         --i;
    //         --j;
    //         // TOREMOVE
    //         if (i<3){
    //             cout<< i<<","<<j<<endl;
    //         }
    //     } else if (bt[i][j] == 1) {
    //         --i;
    //     } else {
    //         --j;
    //     }
    // }

    // int start = j;
    // int end = maxJ;

    // write2DArrayToCSV(kmer_counts,"kmer_counts.csv");
    // write2DArrayToCSV(score,"scores.csv");
    // write2DArrayToCSV(bt,"bt.csv");
    // //////////////////////////////////////////////////////

    return make_pair(start2,end2);
}

int main(int argc, char* argv[]) {

    if (argc != 3 && argc != 6) {
        cerr << "Usage: " << argv[0] << " <fasta A> <fasta B> <match> <mismatch> <indel>\n";
        return 1;
    }
    string fastaA = argv[1];
    string fastaB = argv[2];

    int match = 1;
    int mismatch = -1;
    int gap = -1;
    
    if (argc == 6){
        match = stoi(argv[3]);
        mismatch = stoi(argv[4]);
        gap = stoi(argv[5]);
    }

    bool coordsonly = true; // Output only start and end of the longer string.
    
    string A, B;
    ifstream fileA(fastaA), fileB(fastaB);
    if (!fileA || !fileB) {
        cerr << "Error: Unable to open input files.\n";
        return 1;
    }
    string line;
    getline(fileA, line); // Skip header line
    getline(fileB, line); // Skip header line
    
    while (getline(fileA, line)) {
        A += line;
    }
    while (getline(fileB, line)) {
        B += line;
    }
    string string1 = A;
    string string2 = B;
    // Concatenate strings with $ delimiter
    string concatenated = string1 + "$" + string2;
    
    // Build the suffix array and LCP array
    vector<int> suffix_array;
    vector<int> plcp;
    vector<int> lcp;
    build_SA_LCP(concatenated, suffix_array, plcp, lcp);

    // Find 10-mer hits
    auto [start, end] = find_10mer_hits(concatenated, suffix_array.data(), lcp.data(), string1.size(), string2.size());
    // cout << start << "\n";
    // cout << end << "\n";
}