"""Two classes (single_seq and multi_seq) can be used to study bioinformatics."""

from collections import Counter
from itertools import product
import random

NUC_BASES = ["A", "C", "G", "T"]
REFERENCE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class single_seq:
    """DNA sequnece class. Defalt value: ACGT, No label"""

    def __init__(self, seq="ACGT", label='No label'):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.label = label
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided sequence may not be a correct DNA sequence."

    def __repr__(self):
        return f"[Type]: DNA\n[Sequence]: {self.seq}\n[Length]: {len(self.seq)}\n[Label]: {self.label}"

    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(NUC_BASES).issuperset(self.seq)

    # DNA functions:

    def pattern_count(self, pattern):
        """
        Counts the number of times that a k-mer "pattern" appears as a substring of DNA sequence.

        Parameters:
        ----------
        pattern (str): The k-mer pattern

        Returns
        ----------
        int: The number of times "pattern" appears 
        """
        count = 0
        for i in range(len(self.seq)-len(pattern)+1):
            if self.seq[i:i+len(pattern)] == pattern:
                count += 1
        return count

    def freq_kmers(self, k):
        """
        Finds the most frequent k-mers in DNA sequence.

        Parameters:
        ----------
        k (int): k-mer length

        Returns
        ----------
        list: A list containing all most frequent k-mers
        """
        dict_count = {}
        for i in range(len(self.seq)-k+1):
            if self.seq[i:i+k] in dict_count:
                dict_count[self.seq[i:i+k]] += 1
            else:
                dict_count[self.seq[i:i+k]] = 1
        return [k for k, v in dict_count.items() if v == max(dict_count.values())]

    def reverse_complement(self):
        """
        Finds the transverse complement 

        Returns
        ----------
        str: A string of transverse complement  
        """
        mapping = str.maketrans('ATCG', 'TAGC')
        return self.seq.translate(mapping)[::-1]

    def pattern_index(self, pattern):
        """
        Finds all starting positions where "pattern" appears in the DNA sequence, i.e. pattern matching.

        Parameters:
        ----------
        pattern (str): The k-mer pattern

        Returns
        ----------
        list: A list of integers indicating all starting points     
        """
        index = []
        pl = len(pattern)
        for i in range(len(self.seq)-pl+1):
            if self.seq[i:i+pl] == pattern:
                index.append(i)
        return index

    def clump_find(self, k, L, t):
        """
        Finds all k-mers forming clumps (L,t) in the DNA sequence.

        Parameters:
        ----------
        k (int): k-mer length
        L (int): Moving window size (must be <= len(self.seq))
        t (int): Number of appearance within a given window size

        Returns
        ----------
        set: All distinct k-mers forming (L,t)-clumps 
        """
        match = []
        for _ in range(len(self.seq)-L):
            dict_count = {}
            for i in range(len(self.seq)-k):
                if self.seq[i:i+k] in dict_count:
                    dict_count[self.seq[i:i+k]] += 1
                else:
                    dict_count[self.seq[i:i+k]] = 1
            t_freq = [k for k, v in dict_count.items() if v == t]
            match.extend(t_freq)
        return set(match)

    def min_GCskew(self):
        """
        Find all positions in the DNA sequence where the skew diagram (#G - #C) attains a minimum.

        Returns
        ----------
        set: A list of integers indicating the positions/indices where skew value is minimized
        """
        skew = []
        diff = 0
        for i in self.seq:
            if i == 'G':
                diff += 1
            elif i == 'C':
                diff -= 1
            skew.append(diff)
        minimum = min(skew)
        index = [i for i, x in enumerate(skew) if x == minimum]
        return index

    def approx_pattern_index(self, pattern, dist):
        """
        Finds starting points of all approximate occurrences of a pattern in the DNA sequence, i.e. a generalization of pattern matching.

        Parameters:
        ----------
        pattern (str): The k-mer pattern
        dist (int): The maximum number of mismatches between pattern and a substring of DNA sequence

        Returns
        ----------
        list: A list of integers indicating all starting points 
        """
        index = []
        for loc, i in enumerate(range(len(self.seq)-len(pattern)+1)):
            if HammingDist(self.seq[i:i+len(pattern)], pattern) <= dist:
                index.append(loc)
        return index

    def approx_pattern_count(self, pattern, dist):
        """
        Counts the number of all approximate occurrences of a pattern in the DNA sequence.

        Parameters:
        ----------
        pattern (str): The k-mer pattern
        dist (int): The maximum number of mismatches between pattern and a substring of DNA sequence

        Returns
        ----------
        int: The number of times "approximate" patterns appear
        """
        count = 0
        pl = len(pattern)
        for i in range(len(self.seq)-len(pattern)+1):
            if HammingDist(self.seq[i:i+pl], pattern) <= dist:
                count += 1
        return count

    def freq_approx_kmers(self, k, dist):
        """
        Finds the most frequent k-mers with mismatches in the DNA sequence over all possible k-mers.

        Parameters:
        ----------
        k (int): k-mer length
        dist (int): The maximum number of mismatches between pattern and a substring of DNA sequence

        Returns
        ----------
        list: A list containing all most frequent k-mers with up to "dist" mismatches
        """
        dict_count = {}
        for i in range(len(self.seq)-k+1):
            if self.seq[i:i+k] in dict_count:
                dict_count[self.seq[i:i+k]] += 1
            else:
                dict_count[self.seq[i:i+k]] = 1

        kmers_all = [''.join(base) for base in product(NUC_BASES, repeat=k)]
        ham_count = {}
        for kmer in kmers_all:
            ham_count[kmer] = 0
            for substr in dict_count.keys():
                if HammingDist(kmer, substr) <= dist:
                    ham_count[kmer] += dict_count[substr]
        return [k for k, v in ham_count.items() if v == max(ham_count.values())]

    def freq_approx_rc_kmers(self, k, dist):
        """
        Finds the most frequent k-mers with mismatches and reverse complements in the DNA sequence over all possible k-mers.

        Parameters:
        ----------
        k (int): k-mer length
        dist (int): The maximum number of mismatches between pattern and a substring of DNA sequence

        Returns
        ----------
        list: A list containing all most frequent k-mers with mismatches and reverse complements
        """
        dict_count = {}
        for i in range(len(self.seq)-k+1):
            if self.seq[i:i+k] in dict_count:
                dict_count[self.seq[i:i+k]] += 1
            else:
                dict_count[self.seq[i:i+k]] = 1

        kmers_all = [''.join(base) for base in product(NUC_BASES, repeat=k)]
        ham_count = {}
        for kmer in kmers_all:
            ham_count[kmer] = 0
            for substr in dict_count.keys():
                if HammingDist(kmer, substr) <= dist:
                    ham_count[kmer] += dict_count[substr]
                if HammingDist(RevComplement(kmer), substr) <= dist:
                    ham_count[kmer] += dict_count[substr]
        return [k for k, v in ham_count.items() if v == max(ham_count.values())]

    def profile_kmer(self, k, profile):
        """
        Finds the k-mer that is most likely to be generated by profile among all k-mers in the DNA sequence.

        Parameters:
        ----------
        k (int): k-mer length
        profile (dict): Profile matrix, e.g. {'A': [0.2, 0.2], 'C': [0.4, 0.3], 'G': [0.3, 0.3], 'T': [0.1, 0.2]} for 2-mers.

        Returns
        ----------
        str: A profile-most probable k-mer. 
        NB if there are multiple profile-most probable k-mers, the first instance is returned based on lexicographic order.
        """
        high_prob = 0
        for i in range(len(self.seq)-k+1):
            kmer = self.seq[i:i+k]
            prob = 1
            for w in range(0, len(kmer)):
                prob *= profile[kmer[w]][w]
            if prob > high_prob:
                high_prob = prob
                best_kmer = kmer
        return best_kmer


class multi_seq:
    """A collection of DNA sequneces class. Defalt value: ["ATTTGGC","TGCCTTA","CGGTATC"], No label"""

    def __init__(self, seqs=["ATTTGGC", "TGCCTTA", "CGGTATC"], label='No Label'):
        """Sequence initialization, validation."""
        self.seqs = [seq.upper() for seq in seqs]
        self.label = label
        self.t = len(seqs)  # Number of DNA strings within the list
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided sequence may not be a correct collections of DNA sequences."
        self.valid_length = self.__validlength()
        assert self.valid_length, f"DNA strings have different lengths."

    def __repr__(self):
        return f"[Type]: Collection of DNA strings\n[Sequence]: {self.seqs}\n[Items]: {self.t}\n[Length]: {len(self.seqs[0])}\n[Label]: {self.label}"

    def __validate(self):
        """Check the sequence to make sure it is a valid collection of DNA strings"""
        return all([set(NUC_BASES).issuperset(seq) for seq in self.seqs])

    def __validlength(self):
        return (len(set([len(seq) for seq in self.seqs])) == 1)

    # DNA functions:
    def motif_enumeration(self, k, dist):
        """
        Finds all (k,dist)-motifs in the DNA sequences if the motif appears in every string of DNA with at most d mismatches.

        Parameters:
        ----------
        k (int): motif/k-mer length
        d (int): The maximum number of mismatches between a motif and a substring of DNA sequence

        Returns
        ----------
        list: A list containing all (k,dist)-motifs
        """
        kmers_all = [''.join(base) for base in product(NUC_BASES, repeat=k)]
        kmer_set = []
        for kmer in kmers_all:
            kmatch = 0
            for dnastring in self.seqs:
                for i in range(len(dnastring)-k+1):
                    if HammingDist(kmer, dnastring[i:i+k]) <= dist:
                        kmatch += 1
                        break
            if kmatch == len(self.seqs):
                kmer_set.append(kmer)
        return kmer_set

    def median_string(self, k):
        """
        Finds the medain string that minimizes the Hamming distance between pattern and DNA strings, i.e. min(d(Pattern,DNA)).

        Parameters:
        ----------
        k (int): k-mer length

        Returns
        ----------
        str: A k-mer pattern that minimizes d(Pattern,DNA) among all possible choices of k-mers.
        """
        kmers_all = [''.join(base) for base in product(NUC_BASES, repeat=k)]
        alldist = len(self.seqs)*k
        for kmer in kmers_all:
            matdist = []
            for dnastring in self.seqs:
                rowdist = k
                for i in range(len(dnastring)-k+1):
                    subdist = HammingDist(kmer, dnastring[i:i+k])
                    if subdist < rowdist:
                        rowdist = subdist
                matdist.append(rowdist)
            if sum(matdist) < alldist:
                median = kmer
                alldist = sum(matdist)
        return median

    def greedy_motif_search(self, k, pseudocount=True):
        """
        Finds the profile-most probable motifs/k-mers by greedy search. 
        Specifically, this algorithm starts by forming a motif matrix from arbitrarily selected k-mers in each string from Dna (which in our specific 
        implementation is the first k-mer in each string). It then attempts to improve this initial motif matrix by trying each of the 
        k-mers in Dna1 as the first motif. For a given choice of k-mer Motif1 in Dna1, it builds a profile matrix Profile for this lone k-mer, 
        and sets Motif2 equal to the Profile-most probable k-mer in Dna2. It then iterates by updating Profile as the profile matrix formed 
        from Motif1 and Motif2, and sets Motif3 equal to the Profile-most probable k-mer in Dna3. In general, after finding i − 1 k-mers Motifs 
        in the first i − 1 strings of Dna, greedy_motif_search constructs Profile(Motifs) and selects the Profile-most probable k-mer from Dnai 
        based on this profile matrix. After obtaining a k-mer from each string to obtain a collection Motifs, greedy_motif_search tests to see whether 
        Motifs outscores the current best scoring collection of motifs and then moves Motif1 one symbol over in Dna1, beginning the entire process 
        of generating Motifs again.

        Parameters:
        ----------
        k (int): Motif/k-mer length
        pseudocount (bol, default=True): Apply Laplace's Rule of Succession to form profile matrix

        Returns
        ----------
        list: A collection of motifs/k-mers 
        """
        best_score = 9999
        best_motif = list(map(lambda x: x[0:k], self.seqs))
        for i in range(len(self.seqs[0])-k+1):
            motif = []
            motif.append(self.seqs[0][i:i+k])
            profile = GetLaplaceProfile(
                motif) if pseudocount == True else GetProfile(motif)
            for w in range(1, self.t):
                best_prob = -1
                for j in range(len(self.seqs[1])-k+1):
                    prob = 1
                    for g in range(k):
                        ind = REFERENCE[self.seqs[w][j+g]]
                        prob *= profile[g][ind]
                    if prob > best_prob:
                        best_prob = prob
                        kmer = self.seqs[w][j:j+k]
                motif.append(kmer)
                profile = GetLaplaceProfile(
                    motif) if pseudocount == True else GetProfile(motif)
            score = MedianScore(motif)
            if score < best_score:
                best_score = score
                best_motif = motif
        return best_motif

    def randomized_motif_search(self, k):
        """
        Finds the profile-most probable motifs/k-mers by Monte Carlo algorithm with 1000 iterations. It starts with randomly chosen k-mers 
        in each of t DNA sequences. This algorithm has advantage of being able to find longer motifs.

        Parameters:
        ----------
        k (int): Motif/k-mer length

        Returns
        ----------
        list: A collection of motifs/k-mers 
        """
        max_iter = 1000
        last_motif, last_score = RandomMotifSearch(self.seqs, k, self.t)
        i = 0
        while i < max_iter:
            i += 1
            best_motif, best_score = RandomMotifSearch(self.seqs, k, self.t)
            if best_score < last_score:
                last_motif = best_motif
                last_score = best_score
        return last_motif

    def gibbs_motif_search(self, k, N):
        """
        Finds the profile-most probable motifs/k-mers by Monte Carlo algorithm with 50 iterations. It starts with randomly chosen k-mers 
        in each of t DNA sequences and makes a "random" (based on wewighted profile matrix) rather than a deterministic choice at each iteration.

        Parameters:
        ----------
        k (int): Motif/k-mer length
        N (int): Number of times selecting a profile-randomly generated k-mer at each step 

        Returns
        ----------
        list: A collection of motifs/k-mers 
        """
        max_iter = 30
        last_motif, last_score = GibbsSampler(self.seqs, k, self.t, N)
        i = 0
        while i < max_iter:
            i += 1
            best_motif, best_score = GibbsSampler(self.seqs, k, self.t, N)
            if best_score < last_score:
                last_motif = best_motif
                last_score = best_score
        return last_motif

# Other functions


def HammingDist(seq1, seq2):
    # Calculates the Hamming distance
    count = 0
    for i, v in zip(seq1, seq2):
        if i != v:
            count += 1
    return count


def RevComplement(seq):
    # Returns the reverse complement
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]


def MedianScore(seqs):
    # Calculates the Median score, i.e. the degree of mismatches
    seq_len = len(seqs)
    score = 0
    for i in range(len(seqs[0])):
        nucs = list(map(lambda x: x[i], seqs))
        counter = Counter(nucs)
        score += seq_len - max(counter.values())
    return score


def GetProfile(seq):
    # Caculates the profile matrix based on a collection of strings
    profile = []
    seq_len = len(seq)
    for i in range(len(seq[0])):
        holder = list(map(lambda x: x[i], seq))
        ratio = []
        for b in NUC_BASES:
            ratio.append(holder.count(b)/seq_len)
        profile.append(ratio)
    return profile


def GetLaplaceProfile(motif):
    # Caculates the profile matrix based on a collection of strings using pseudocounts
    profile = []
    seq_len = len(motif)
    for i in range(len(motif[0])):
        holder = list(map(lambda x: x[i], motif))
        ratio = []
        for b in NUC_BASES:
            ratio.append((holder.count(b)+1)/(seq_len+4))
        profile.append(ratio)
    return profile


def PatternProfileSeq(seq, k, profile):
    # Deterministically choose a k-mer based on the weights in profile matrix
    best_prob = 0
    for j in range(len(seq)-k+1):
        prob = 1
        for g in range(k):
            ind = REFERENCE[seq[j+g]]
            prob *= profile[g][ind]
        if prob > best_prob:
            best_prob = prob
            kmer = seq[j:j+k]
    return kmer


def RandomMotifSearch(seqs, k, t):
    # Randomized Motif Search
    seq_len = len(seqs[0])
    motif = []
    for seq in seqs:
        rand_start = random.randrange(seq_len-k+1)
        motif.append(seq[rand_start:rand_start+k])

    best_motif = motif
    best_score = MedianScore(motif)
    while True:
        profile = GetLaplaceProfile(motif)
        motif = []
        for f in range(t):
            motif_iter = PatternProfileSeq(seqs[f], k, profile)
            motif.append(motif_iter)

        score = MedianScore(motif)
        if score < best_score:
            best_motif = motif
            best_score = score
        else:
            return best_motif, best_score


def RandProfileKmer(seq, k, profile):
    # Randomly choose a k-mer based on the weights in profile matrix
    N = len(seq)-k+1
    weights = []
    for j in range(N):
        prob = 1
        for g in range(k):
            ind = REFERENCE[seq[j+g]]
            prob *= profile[g][ind]
        weights.append(prob)
    irand = random.choices(list(range(N)), weights=weights, k=1)[0]
    return seq[irand:irand+k]


def GibbsSampler(seqs, k, t, N):
    # Gibbs Sampler
    seq_len = len(seqs[0])
    motif = []
    for seq in seqs:
        rand_start = random.randrange(seq_len-k+1)
        motif.append(seq[rand_start:rand_start+k])
    best_motif = motif
    best_score = MedianScore(motif)

    for _ in range(N):
        i_rand = random.randrange(t)
        temp = motif.copy()
        temp.pop(i_rand)
        profile = GetLaplaceProfile(temp)
        i_motif = RandProfileKmer(seqs[i_rand], k, profile)
        motif[i_rand] = i_motif
        score = MedianScore(motif)
        if score < best_score:
            best_motif = motif
            best_score = score
    return best_motif, best_score
