from dna_seq import single_seq, multi_seq

# Examples of single_seq
test_dna = single_seq(
    "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA")
print(test_dna)
print(test_dna.pattern_count("AGA"))
print(test_dna.freq_kmers(4))
print(test_dna.pattern_index("CGA"))
print(test_dna.clump_find(5, 50, 4))
print(test_dna.min_GCskew())
print(test_dna.approx_pattern_index("ATTCTGGA", 3))
print(test_dna.approx_pattern_count("GAGG", 2))
print(test_dna.freq_approx_kmers(4, 1))
profile = {
    "A": [0.2, 0.2, 0.3, 0.2, 0.3],
    "C": [0.4, 0.3, 0.1, 0.5, 0.1],
    "G": [0.3, 0.3, 0.5, 0.2, 0.4],
    "T": [0.1, 0.2, 0.1, 0.1, 0.2]}
print(test_dna.profile_kmer(5, profile))

print("-"*20)
# # Examples of multi_seq
test_dna2 = multi_seq(["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                       "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", 
                       "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"])
print(test_dna2)
print(test_dna2.motif_enumeration(5, 1))
print(test_dna2.median_string(6))
print(test_dna2.greedy_motif_search(3, False))
print(test_dna2.randomized_motif_search(8))
print(test_dna2.gibbs_motif_search(8, 100))
