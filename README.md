# gene_sequencing
Code implements two versions of a dynamic programming algorithm for computing the minimal cost of aligning gene sequences and extracting an optimal alignment. 
Dynamic Programming is mainly an optimization over plain recursion. Wherever we see a recursive solution that has repeated calls for same inputs, 
we can optimize it using Dynamic Programming. The idea is to simply store the results of subproblems, so that we do not have 
to re-compute them when needed later. This simple optimization reduces time complexities from exponential to polynomial.
