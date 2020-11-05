"""Two DNA classes (dna_fragments and dna_pairs) can be used to reconstruct a DNA seqeunce from parts."""

from bio_table import NUC_BASES
from collections import Counter
import random

class dna_fragments:
    """A collection of DNA strings class. Default value: ['ACCGA','CCGAA','CGAAG'], No label"""

    def __init__(self, seqs=['ACCGA','CCGAA','CGAAG'], label='No Label'):
        """Sequence initialization, validation."""
        self.seqs = [seq.upper() for seq in seqs]
        self.label = label
        self.t = len(seqs)  # Number of DNA strings within the list
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided sequence may not be a correct collections of DNA strings."
        self.valid_length = self.__validlength()
        assert self.valid_length, f"DNA strings have different lengths."
    
    def __validate(self):
        """Check the sequence to make sure it is a valid collection of DNA strings"""
        return all([set(NUC_BASES).issuperset(seq) for seq in self.seqs])

    def __validlength(self):
        return (len(set([len(seq) for seq in self.seqs])) == 1)

    # DNA functions:
    def genome_path(self):
        """
        Reconstructs the genome from its genome path where the fragments are in the right order from left to right.

        Returns
        ----------
        str: The reconstructed genome
        """
        return GenomePath(self.seqs)

    def overlap_graph(self):
        """
        Constructs the overlap graph of a collection of k-mers.

        Returns
        ----------
        dict: An overlap (direced) graph, in the form of an adjacency list
        """
        uni_seqs = set(self.seqs)
        Overlap = {}
        for seq in uni_seqs:
            pattern = seq[1:]
            for temp in uni_seqs:
                if pattern == temp[:-1]:
                    if seq not in Overlap:
                        Overlap[seq] = [temp]
                    else:
                        Overlap[seq].append(temp)
        return Overlap

    def debruijn_graph(self):
        """
        Constructs the De bruijn graph of a collection of k-mers.

        Returns
        ----------
        dict: A De bruijn graph, in the form of an adjacency list
        """
        return debruijnSeq(self.seqs)

    def string_construct(self):
        """
        Reconstructs the genome from a list of k-mers patterns.

        Returns
        ----------
        str: The reconstructed genome
        """
        db = debruijnSeq(self.seqs)
        path = EulerianPath(db)
        text = GenomePath(path)
        return text

    def find_contigs(self):
        """
        Finds all contigs (long, contiguous segments of the genome) .

        Returns
        ----------
        list: All contigs in DeBruijn(Patterns)
        """
        graph = debruijnSeq(self.seqs)
        a = MaximalNonBranchingPaths(graph)
        return list(map(GenomePath,a))
    

class dna_pairs:
    """Paired Reads, i.e. (k,d)-mer class. Default value: ['GACC|GCGC','ACCG|CGCC','CCGA|GCCG','CGAG|CCGG','GAGC|CGGA'], No label"""

    def __init__(self, seqs=['GACC|GCGC','ACCG|CGCC','CCGA|GCCG','CGAG|CCGG','GAGC|CGGA'], label='No Label'):
        """Sequence initialization, validation."""
        self.seqs = [seq.upper() for seq in seqs]
        self.label = label
        self.t = len(seqs)  # Number of pairs within the list
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided sequence may not be a correct collections of DNA reads or in the right format."
        self.valid_length = self.__validlength()
        assert self.valid_length, f"DNA strings have different lengths."
    
    def __validate(self):
        """Check the sequence to make sure it is a valid collection of DNA strings"""
        return all([set(NUC_BASES+['|']).issuperset(seq) for seq in self.seqs])

    def __validlength(self):
        return (len(set([len(seq) for seq in self.seqs])) == 1)

    # DNA functions:
    def string_construct(self,k,d):
        """
        Reconstructs the genome from a collection of paired (k,d)-mers PairedReads.

        Parameters:
        ----------
        k (int): k-mer length
        d (int): gap length between reads

        Returns
        ----------
        str: The reconstructed genome
        """
        pairs = {}
        for seq in self.seqs:
            pair = seq.split('|')
            pre = '|'.join(map(lambda x: x[:-1],pair))
            post = ['|'.join(map(lambda x: x[1:],pair))]
            if pre not in pairs:
                pairs[pre] = post
            else:
                pairs[pre].append(post)
        gapped = EulerianPath(pairs)
        # for node in gapped_pattern:
        first = GenomePath(list(map(lambda x: x.split('|')[0],gapped)))
        second = GenomePath(list(map(lambda x: x.split('|')[1],gapped)))
        if first[k+d:] == second[:-(k+d)]:
            pair_read = first+second[-(k+d):] 
            return pair_read
        else: 
            return "THERE IS NO STRING SPELLED BY THE GAPPED PATTERNS"


# Other functions

def GenomePath(seqs):
    # Reconstructs the genome from its genome path, which is in the right order
    lastnuc = ''.join(map(lambda x: x[-1],seqs))
    return seqs[0]+lastnuc[1:]

def debruijnSeq(seqs):
    # Constructs the De bruijn graph of a collection of k-mers
    Adjacency = {}
    for seq in seqs:
        pattern = seq[:-1]
        temp = seq[1:] 
        if pattern not in Adjacency:
            Adjacency[pattern] = [temp]
        else:
            Adjacency[pattern].append(temp)
    return Adjacency

def EulerianCycle(graph):
    # Creates an Eulerian cycle from an adjacency list 
    stack = []
    circuit = []
    location = random.choice(list(graph.keys())) #random start
    while True:
        if location in graph.keys():
            nextloc = random.choice(graph[location])
            stack.append(location)
            if len(graph[location]) > 1:
                graph[location].remove(nextloc)
            else:
                del graph[location]
            location = nextloc
        else:
            circuit.append(location)
            location = stack.pop()
        if graph == {}:
            return stack+[location]+circuit[::-1]

def EulerianPath(graph):
    # Creates an Eulerian path from an adjacency list 
    import random
    stack = []
    circuit = []
    location = SubstractDegree(graph)
    while True:
        if location in graph.keys():
            nextloc = random.choice(graph[location])
            stack.append(location)
            if len(graph[location]) > 1:
                graph[location].remove(nextloc)
            else:
                del graph[location]
            location = nextloc
        else:
            circuit.append(location)
            location = stack.pop()
        if graph == {}:
            return stack+[location]+circuit[::-1]

def SubstractDegree(graph):
    # Calculate the difference between in-degree and out-degree of a node
    graph_in= []
    graph_out = {}
    for key, val in graph.items():
        graph_in.extend(val)
        graph_out[key] = len(val)
    graph_in = dict(Counter(graph_in))
    subdict = {key: graph_out[key] - graph_in.get(key, 0) for key in graph_out}
    for key,i in subdict.items():
         if i==1: location = key
    return location

def DegreeCalc(graph):
    # Calculate in-degree and out-degree of a node
    all_values = []
    for i in graph.values():
        all_values.extend(i)
    keys = list(graph.keys())
    all_keys= set(all_values+keys) 
    graph_in = []
    for val in graph.values():
        graph_in.extend(val)
    graph_in = dict(Counter(graph_in))
    degree = {}
    for item in all_keys:
        indegree = graph_in[item] if item in graph_in else 0
        outdegree = len(graph[item]) if item in graph else 0
        degree[item]= [indegree,outdegree]
    return degree
    
def MaximalNonBranchingPaths(graph):
    # Finds all maximal non-branching paths in a De bruijn graph
    degree = DegreeCalc(graph)
    paths = []
    for key, val in degree.items():
        if val != [1,1]:
            if val[1] > 0:
                for conn in graph[key]:
                    path = [key,conn]
                    nextloc = conn                   
                    din,dout = degree[nextloc]                    
                    while (din==1) and (dout ==1):
                        temp = graph[nextloc][0]
                        del graph[nextloc]
                        path.append(temp)
                        nextloc = temp
                        din,dout = degree[nextloc]
                    else:
                        paths.append(path)
                del graph[key]
    if graph !={}:
        paths.append(EulerianCycle(graph))          
    return paths

# Print graph function
def print_graph(graph):
    for key in graph.keys():
        print(key,'-> ',end='')
        print(*graph[key],sep=',')