from dna_seq import single_seq, multi_seq
from dna_assembly import print_graph

print('Examples of single_seq class:')

test_dna = single_seq(
    "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA")

# Returns the object representation
print(test_dna)

# Counts the number of times that a k-mer "pattern" appears as a substring
print(test_dna.pattern_count("AGA"))
# Output: 7

# Finds the most frequent k-mers in DNA sequence.
print(test_dna.freq_kmers(4))
# Output: ['CGAC', 'GACA', 'GAAG', 'AAGA']

# Finds the transverse complement
print(test_dna.pattern_index("CGA"))
# Output: [6, 22, 38, 43]

# Finds all k-mers forming clumps (L,t) in the DNA sequence
print(test_dna.clump_find(5, 50, 4))
# Output: {'GAAGA', 'CGACA'}

# Find all positions in the DNA sequence where the skew diagram (#G - #C) attains a minimum
print(test_dna.min_GCskew())
# Output: [0, 6, 9, 10]

# Finds starting points of all approximate occurrences of a pattern in the DNA sequence, i.e. a generalization of pattern matching
print(test_dna.approx_pattern_index("ATTCTGGA", 3))
# Output: [26, 33]

#  Counts the number of all approximate occurrences of a pattern in the DNA sequence
print(test_dna.approx_pattern_count("GAGG", 2))
# Output: 28

# Finds the most frequent k-mers with mismatches in the DNA sequence over all possible k-mers
print(test_dna.freq_approx_kmers(4, 1))
# Output: ['GAGA']

profile = {
    "A": [0.2, 0.2, 0.3, 0.2, 0.3],
    "C": [0.4, 0.3, 0.1, 0.5, 0.1],
    "G": [0.3, 0.3, 0.5, 0.2, 0.4],
    "T": [0.1, 0.2, 0.1, 0.1, 0.2]}
print(test_dna.profile_kmer(5, profile))
# Output: CGACA

# Translate (& transcribe) the DNA sequence into an amino acid string
test_dna = single_seq("ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA")
print(test_dna.translation())
# Output: MAMAPRTEINSTRING

# Find substrings of then DNA sequence encoding a given amino acid sequence
print(test_dna.find_codon_string('MA'))
# Output: ['ATGGCC', 'ATGGCC', 'GGCCAT']

# Generates the k-mer composition of the DNA sequence
test_dna = single_seq("CAATCCAAC")
print(test_dna.construct_kmer(5))
# Output: {'CAATC', 'AATCC', 'ATCCA', 'TCCAA', 'CCAAC'}

# Constructs a deBruijn graph of the DNA sequence
test_dna = single_seq("AAGATTCTCTAAGA")
print_graph(test_dna.debrujin_graph((4)))
# Output: 
# AAG -> AGA,AGA
# AGA -> GAT
# GAT -> ATT
# ATT -> TTC
# TTC -> TCT
# TCT -> CTC,CTA
# CTC -> TCT
# CTA -> TAA
# TAA -> AAG

print("-"*20)

print('Examples of multi_seq class:')

test_dna2 = multi_seq(["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                       "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", 
                       "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"])
# Returns the object representation                      
print(test_dna2)

# Finds all (k,dist)-motifs in the DNA sequences if the motif appears in every string of DNA with at most d mismatches
print(test_dna2.motif_enumeration(5, 1))
# Output: ['AGGTG', 'AGTGT', 'AGTTC', 'GTGTA', 'GTTCC', 'TCCAG']

# Finds the medain string that minimizes the Hamming distance between pattern and DNA strings, i.e. min(d(Pattern,DNA))
print(test_dna2.median_string(6))
# Output: AGGTGC

# Finds the profile-most probable motifs/k-mers by greedy search
print(test_dna2.greedy_motif_search(3, False))
# Output: ['CAG', 'CAG', 'CAG', 'CAG', 'CAG']

# Finds the profile-most probable motifs/k-mers by Monte Carlo algorithm with 1000 iterations
print(test_dna2.randomized_motif_search(8))
# Output: ['AACGGCCA', 'AAGTGCCA', 'TAGTACCG', 'AAGTTTCA', 'ACGTGCAA']

# Finds the profile-most probable motifs/k-mers by Monte Carlo algorithm with 50 iterations
print(test_dna2.gibbs_motif_search(8, 100))
# Output: ['TAAACGGC', 'TGTAAGTG', 'TACAGGCG', 'TGCACGTC', 'TCCACGTG']

print("-"*20)