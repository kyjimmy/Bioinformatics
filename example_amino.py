from amino_acid import peptide, mass_spectrum

print('Examples of peptide class:')

# Calculates the integer mass of the peptide
test_dna = peptide('LEQN')
print(test_dna.mass()) 
# Output: 484

# Generates the theoretical spectrum of a cyclic peptide
print(test_dna.cyclospectrum()) 
# Output: [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]

# Counts the number of masses shared between Cyclospectrum(peptide) and a given spectrum
test_dna = peptide('NQEL')
print(test_dna.score([0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484])) 
# Output: 11

print('-'*20)

print('Examples of mass_spectrum class:')

# Finds peptides that are most consistent with the spectrum 
test_dna = mass_spectrum([0, 113, 128, 186, 241, 299, 314, 427])
for i in test_dna.cyclopeptide_sequencing(): print(*i,sep='-')
# Output: 113-128-186 113-186-128 128-113-186 128-186-113 186-113-128 186-128-113

# Finds peptides that are most consistent with the spectrum using Leaderboard Cyclopeptide Sequencing
test_dna = mass_spectrum([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460])
print(*test_dna.leaderboard_cyclopeptide_sequencing(10)[0],sep='-') 
# Output: 71-147-113-129

# Finds peptides that are most consistent with the spectrum using Spectral Convolution and Leaderboard Cyclopeptide Sequencing
test_dna = mass_spectrum([57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493])
print(*test_dna.convo_cyclopeptide_sequencing(20,60)[0],sep='-') 
# Output: 57-129-99-71-80-57

print('-'*20)