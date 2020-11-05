from dna_assembly import dna_fragments, dna_pairs, print_graph

print('Examples of dna_fragments class:')

# Reconstructs the genome from its genome path
test_dna = dna_fragments(['ACCGA','CCGAA','CGAAG','GAAGC','AAGCT'])
print(test_dna.genome_path())
# Output: ACCGAAGCT

# Constructs the overlap graph
test_dna = dna_fragments(['ATGCG','GCATG','CATGC','AGGCA','GGCAT','GGCAC'])
print_graph(test_dna.overlap_graph())
# Output: 
# AGGCA -> GGCAT,GGCAC
# GCATG -> CATGC
# GGCAT -> GCATG
# CATGC -> ATGCG

# Constructs the De bruijn graph
test_dna = dna_fragments(['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG'])
print_graph(test_dna.debruijn_graph())
# Output: 
# GAG -> AGG
# CAG -> AGG,AGG
# GGG -> GGG,GGA
# AGG -> GGG
# GGA -> GAG

# Reconstructs the genome
test_dna = dna_fragments(['CTTA','ACCA','TACC','GGCT','GCTT','TTAC'])
print(test_dna.string_construct())
# Output: GGCTTACCA

# Finds all contigs
test_dna = dna_fragments(['ATG','ATG','TGT','TGG','CAT','GGA','GAT','AGA'])
print(test_dna.find_contigs())
# Output: ['CAT', 'AGA', 'ATG', 'ATG', 'GAT', 'TGT', 'TGGA']

print('-'*20)

print('Examples of dna_pairs class:')

# Reconstructs the genome from a collection of paired (k,d)-mers PairedReads
test_dna = dna_pairs(['GACC|GCGC','ACCG|CGCC','CCGA|GCCG','CGAG|CCGG','GAGC|CGGA'])
print(test_dna.string_construct(4,2))
# Output: GACCGAGCGCCGGA

print('-'*20)