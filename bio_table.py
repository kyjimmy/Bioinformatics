NUC_BASES = ['A', 'C', 'G', 'T']

AMINO_BASES = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

DNA_CODONS = {
    # 'M': START CODON, '_': STOP CODON
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A", 
    "TGC": "C", "TGT": "C", 
    "GAC": "D", "GAT": "D",
    "GAA": "E", "GAG": "E",
    "TTC": "F", "TTT": "F", 
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", 
    "CAC": "H", "CAT": "H", 
    "ATA": "I", "ATC": "I", "ATT": "I",
    "AAA": "K", "AAG": "K",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", "TTA": "L", "TTG": "L", 
    "ATG": "M",
    "AAC": "N", "AAT": "N", 
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P", 
    "CAA": "Q", "CAG": "Q",
    "AGA": "R", "AGG": "R", "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", 
    "AGC": "S", "AGT": "S", "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T", 
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"}
    
AMINO_MASS = {
    "G": 57, "A": 71, "S": 87, "P": 97, "V": 99, "T": 101, "C": 103, "I": 113,
    "L": 113, "N": 114, "D": 115, "K": 128, "Q": 128, "E": 129, "M": 131,
    "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186}

PROTEINOGENIC_MASS = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]