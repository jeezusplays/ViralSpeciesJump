#Dijkstra

from os import path
from heapq import heappush, heappop
from scipy.spatial import distance 
import time
import pandas as pd
#import networkx as nx
#import matplotlib as plt

wild = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"
len_wild = len(wild) #get str length of wildtype

num_snp_col = 'num_snp' #number of substitutions
seq_data_col = 'aa_seq' #converted to amino acid sequence

codons = { #codon dictionary
    'F': ["TTT","TTC"],
    'L': ["TTA","TTG", "CTT","CTC", "CTA", "CTG"],
    'I': ["ATT","ATC", "ATA"],
    'M': ["ATG"],
    'V': ["GTT","GTC", "GTA", "GTG"],
    'S': ["TCT","TCC", "TCA", "TCG", "AGT","AGC"],
    'P': ["CCT","CCC", "CCA", "CCG"],
    'T': ["ACT","ACC", "ACA", "ACG"],
    'A': ["GCT","GCC", "GCA", "GCG"],
    'Y': ["TAT","TAC"],
    "O": ["TAA", "TAG", "TGA"], #STOP codon
    'H': ["CAT","CAC"],
    'Q': ["CAA","CAG"],
    'N': ["AAT","AAC"],
    'K': ["AAA","AAG"],
    'D': ["GAT","GAC"],
    'E': ["GAA","GAG"],
    'C': ["TGT","TGC"],
    'W': ["TGG"],
    'R': ["CGT","CGC", "CGA", "CGG", "AGA", "AGG"],
    'G': ["GGT","GGC", "GGA", "GGG"]
}

# START =========== process excel data =============

def process_sequences(): #feature to process original data
    df = pd.read_csv(r"../Data/CSV_PR_original.csv") #reads csv
    df = df[~df.NASeq.str.contains('~')] #remove rows with ~ in sequences
    mask = df.NASeq.map(len) == len_wild #retrieve rows where sequences are same length as wildtype
    df = df[mask].reset_index() #new df, index reset
    return df.NASeq

def process_data(raw = False): 
    seq_data = process_sequences()

    if raw == False: #if raw sequence is false
        wild_aa = ''.join(aa_gen(wild)) #create wildtype amino acid
        seq_data = pd.df(aa_gen(seq_data)) #create amino sequences from data sequences
        seq_data.columns = [seq_data_col] 
        mask = seq_data[seq_data_col].map(len) == len(wild_aa) #both aa sequence same len if not remove
        seq_data = seq_data[mask].reset_index()    

    c = [] #empty list
    for seq in seq_data[seq_data_col]:
        c.append(distance.levenshtein(wild_aa, seq)) #calculate levenshtein distance and add to end of empty list c

    c = pd.df(c) #make list c a df
    
    df = pd.concat([c[0], seq_data], axis=1) #create new df by combining hamming distance and aa sequence
    df.columns = [num_snp_col, 'index', seq_data_col] 
    df = df[['index', num_snp_col, seq_data_col]] #re order columns
    df.drop_duplicates([seq_data_col], keep="first", inplace=True) #drop duplicates
    return df

def process_SNP_probabilities(): #feature to process the prob of change
    df = pd.read_csv(r"../Data/probabilties.csv")
    df.set_index('AA', inplace=True)
    df = df.round(3) #round df to 3 decimal places
    return df

def get_aa_data(filename): #retrieve aa sequences
    if path.exists("../Data/" + filename):
        print("Found data file")
        data = pd.read_csv(r"../Data/" + filename)#.head(max)
    else:
        data = process_data(False)
        data.to_csv(r"../Data/" + filename, index=False)
    return data

def get_probs(): #retrieve probability
    if path.exists("../Data/probs.csv"):
        print("Found probs file")
        probs = pd.read_csv(r"../Data/probs.csv").set_index('AA')
    else:
        probs = process_SNP_probabilities()
        probs.to_csv(r"../Data/probs.csv")
    return probs

def print_groups(groups):
    for name_of_the_group, group in groups:
        print ("Group = " + str(int(name_of_the_group)))
        print (group)

# END ============ process excel data ===============

# START ===== nucleotide to amino acid sequence =====

def aa_gen(sequence): #amino acid generation
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        aa = codon_aa(codon.upper())
        if aa == None:
            print("aa_gen(): codon error " + codon)
            continue
        else:
            yield aa

def codon_aa(codon): #retrieve number of substitutions = keys
    for key, value in codons.items():
        if(codon in value):
            return key

# END ===== nucleotide to amino acid sequence =====

# START =========== child processing ==============

def get_SNPs(s1, s2):
    indices = [i for i in range(len(s1)) if s1[i] != s2[i]] #list comprehension for each character check 
    snps = []
    for i in indices:
        snps.append({'index':i, 'codon1':s1[i], 'codon2':s2[i]}) #dictionary of the mutations from s1 to s2 with the index
    return snps

def get_seq_from_group(groups, group_num, index_in_group):
    if group_num in groups.groups.keys():
        group = get_group(groups, group_num)
        return group.iloc[index_in_group, 2]
    else:
        print("Group not found")
        return None

def get_group(groups, group_num):
    if group_num in groups.groups.keys():
        return groups.get_group(group_num)
    else:
        print("Group not found")
        return None
    
# END =========== child processing ==============

# START ========== Dijkstra's Algorithm =================

#def calc_h(node1, node2): #heuristic function
    #n = 0 #number of 0 mutations
    #p = 1
    #muts = [[x,y] for x,y in zip(node1.sequence, node2.sequence) if x != y]
    #for m in muts:
        #mp = probs.loc[m[0], m[1]]
        #if mp == 0: #mutation = 0
            #n = n + 1
        #else:
            #p = p * mp
    #return (1/p) + n

class PriorityQueue:
  def __init__(self):
    self.pq = []

  def print(self):
      for n in self.pq:
          print(n.sequence)

  def add(self, item, priority):
    heappush(self.pq, (priority, item))

  def poll(self):
    return heappop(self.pq)

  def peek(self):
    return self.pq[0]

  def remove(self, item):
    value = self.pq.remove(item)
    heapify(self.pq)
    return value is not None

  def __len__(self):
    return len(self.pq)

class Node:
    def __init__(self, sequence, depth=0, parent=None):
        self.sequence = sequence
        self.hash = hash(self.sequence)
        self.parent = parent
        self.depth = depth
        #self.hcost = 0
        self.gcost = self.g()

    #def levenshtein(self):
        #seq_len = len(self.sequence)
        #return len([i for i in range(seq_len) if self.sequence[i] != self.parent.sequence[i]])

    def get_neighbours(self):
        for key in list(keys)[self.depth+1:]:
            bin = bins.get_group(key)#.reset_index()
            neighbours = [Node(seq, self.depth+1, self) for seq in bin[seq_data_col]]
            return neighbours

    def f(self):
        return  self.gcost #+ self.hcost

    def g(self):
        gcost = 0
        if self.parent != None: #no parent, gcost = 0
            gcost = self.parent.gcost #+ calc_h(self.parent, self) #got parent, take the parent gcost and calc cost via prob
        else:
            gcost = 0
        return gcost

    def __cmp__(self, other): #compare function
        if not isinstance(other, Node):
            return NotImplemented
        return self.sequence == other.sequence

    def __eq__(self, other):
        return self.__cmp__(other)

    def __hash__(self):
        return self.hash #hash(self.sequence)

    def __lt__(self, other): #less than function
        return self.f() < other.f()

class Solver:
    def rebuild_path(self, end):
        path = [end]
        node = end.parent
        while True:
          path.append(node)
          if node.parent == None:
              break
          node = node.parent
        return path

    def solve(self, start, goal):
        max_n = 20 #limit neighbours
        open = PriorityQueue()
        open.add(start, 0)
        closed = set()
        while open:
            current = open.poll()

            if current.sequence == goal.sequence: ###tuple object has no attribute 'sequence'
                print("----PATH----")
                path = self.rebuild_path(current)
                for node in reversed(path):
                    print(node.sequence)
                break
            neighbours = current.get_neighbours()
            
            n = 0
            for node in neighbours:
                #node.hcost = calc_h(node, goal)
                if node not in closed:
                    open.add(node)
                n = n + 1
                if n >= max_n:
                    break
            closed.add(current)

# END ================= A* ====================

data = get_aa_data("data_aa.csv") 
#data = get_raw_data().reset_index()

data.drop_duplicates([seq_data_col], keep="first", inplace=True)
probs = get_probs().transpose() #row become column, column become row
bins = data.groupby(num_snp_col)

print_groups(bins)
keys = bins.groups.keys()
#print(keys)

wild = ''.join(aa_gen(wild))
start = Node(wild, 0, None)

start = Node(get_seq_from_group(bins, 0, 0), 0, None)
goal = Node(get_seq_from_group(bins, 7, 0), 7, None)

start_time = time.process_time()

# dict_keys([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,31])
depth = 3
index = 0
parent = None
sequence = get_seq_from_group(bins, depth, index)
goal = Node(sequence, depth, parent)

print("-----------")
print("Start: " + start.sequence)
print("End: " + goal.sequence)

if goal:
    astar = Solver()
    astar.solve(start, goal)

print ("Execution Time: " + str(time.process_time() - start_time))