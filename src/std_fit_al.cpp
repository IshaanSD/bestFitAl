#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

int max3(int a, int b, int c) {
    return max(max(a, b), c);
}

// Function to compute the fitting alignment score matrix and backtrack
pair<vector<vector<int>>, vector<vector<int>>> fittingAlignment(const string& A, const string& B, int match, int mismatch, int gap) {
    int m = A.size();
    int n = B.size();

    // Initialize the score matrix and the backtrack matrix
    vector<vector<int>> score(m + 1, vector<int>(n + 1, 0));
    vector<vector<int>> bt(m + 1, vector<int>(n + 1, 0));

    // Fill the score matrix
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int matchMismatch = score[i - 1][j - 1] + (A[i - 1] == B[j - 1] ? match : mismatch);
            int deleteFromA = score[i - 1][j] + gap;
            int insertToA = score[i][j - 1] + gap;
            if (score[i][j - 1] == 0) insertToA = max(0,insertToA);
            score[i][j] = max3(matchMismatch, deleteFromA, insertToA);

            // Track the direction from which the score was derived
            if (score[i][j] == matchMismatch) {
                bt[i][j] = 0; // Diagonal
            } else if (score[i][j] == insertToA) {
                bt[i][j] = 1; // left
            } else {
                bt[i][j] = 2; // up
            }
        }
    }

    return {score, bt};
}

pair<int, int> backtrackCoords(const string& A, const string& B, const vector<vector<int>>& score, const vector<vector<int>>& bt) {
    int m = A.size();
    int n = B.size();

    int maxScore = score[m][0];
    int maxJ = 0;
    for (int j = 1; j <= n; ++j) {
        if (score[m][j] > maxScore) {
            maxScore = score[m][j];
            maxJ = j;
        }
    }

    int i = m, j = maxJ;
    while (i > 0 && j > 0) {
        if (bt[i][j] == 0) {
            --i;
            --j;
        } else if (bt[i][j] == 1) {
            --i;
        } else {
            --j;
        }
    }

    int start = j;
    int end = maxJ;

    return {start, end};
}

// Function to perform the backtrack and get the aligned sequences
pair<string, string> backtrack(const string& A, const string& B, const vector<vector<int>>& score, const vector<vector<int>>& bt) {
    int m = A.size();
    int n = B.size();

    // Find the position of the maximum score in the last row
    int maxScore = score[m][0];
    int maxJ = 0;
    for (int j = 1; j <= n; ++j) {
        if (score[m][j] > maxScore) {
            maxScore = score[m][j];
            maxJ = j;
        }
    }

    // Backtrack to get the aligned sequences
    string alignedA, alignedB;
    int i = m, j = maxJ;
    while (i > 0 && j > 0) {
        if (bt[i][j] == 0) {
            alignedA = A[i - 1] + alignedA;
            alignedB = B[j - 1] + alignedB;
            --i;
            --j;
        } else if (bt[i][j] == 2) {
            alignedA = A[i - 1] + alignedA;
            alignedB = '-' + alignedB;
            --i;
        } else {
            alignedA = '-' + alignedA;
            alignedB = B[j - 1] + alignedB;
            --j;
        }
    }

    return {alignedA, alignedB};
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
    getline(fileA, A); // Skip header line
    getline(fileB, B); // Skip header line
    string line;
    while (getline(fileA, line)) {
        A += line;
    }
    while (getline(fileB, line)) {
        B += line;
    }

    auto [score, bt] = fittingAlignment(A, B, match, mismatch, gap);
    if (!coordsonly){
        auto [alignedA, alignedB] = backtrack(A, B, score, bt);
        ////////////////////////////////////////////////////////////////
        // Output matrices (debug mode)
        // cout << "Alignment scores matrix:\n";
        // cout << "\t\t";
        // for (char b : B) {
        //     cout << b << "\t";
        // }
        // cout << "\n\t";
        // for (size_t j = 0; j <= B.size(); ++j) {
        //     cout << score[0][j] << "\t";
        // }
        // cout << "\n";
        
        // for (size_t i = 1; i <= A.size(); ++i) {
        //     cout << A[i - 1] << "\t";
        //     for (size_t j = 0; j <= B.size(); ++j) {
        //         cout << score[i][j] << "\t";
        //     }
        //     cout << "\n";
        // }
        
        // cout << "Backtrack matrix:\n";
        // cout << "\t\t";
        // for (char b : B) {
        //     cout << b << "\t";
        // }
        // cout << "\n\t";
        // for (size_t j = 0; j <= B.size(); ++j) {
        //     cout << bt[0][j] << "\t";
        // }
        // cout << "\n";
        
        // for (size_t i = 1; i <= A.size(); ++i) {
        //     cout << A[i - 1] << "\t";
        //     for (size_t j = 0; j <= B.size(); ++j) {
        //         cout << bt[i][j] << "\t";
        //     }
        //     cout << "\n";
        // }
        ////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////
        // Main output - aligned strings
        cout << alignedA << "\n";
        cout << alignedB << "\n";
        ////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////
        // Main output - aligned strings
        cout << alignedA << "\n";
        cout << alignedB << "\n";
        ////////////////////////////////////////////////////////////////
    }
    else{
        auto [start, end] = backtrackCoords(A, B, score, bt);
        cout << start << "\n";
        cout << end << "\n";
    }



    return 0;
}