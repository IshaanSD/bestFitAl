# BestFitAl
## Approximate overlapping alignment
For two sequences P (longer) and Q, we find all shared substrings longer than 10 bp between the two using suffix and longest common prefix arrays. We count the number of shared substrings (x,y) starting from position x in P and position y in Q lying on the same diagonal (deonted by x-y). We aggregate diagonals using windows of size B=100bp (deonted by floor((x-y)/B)) and find the diagonal window with highest number of shared substrings. The start and end of this window in the two sequences (thus cropping one or both of the sequences) are used as new boundaries for 2D block-rectangle in the scaled blockplot.

## Standard Fitting Alignment
Implemented using the description available at https://rosalind.info/glossary/fitting-alignment/
