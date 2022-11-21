#A* Algorithm, Base Notation Probability (A, C, G, T)
from os import path
from heapq import heappush, heappop
from scipy.spatial import distance 
import time
import pandas as pd

wild = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"

sub_col = 'num_snp' #number of substitutions
seq_col = 'aa_seq' #converted to amino acid sequence

codons = { #codon dictionary
    'C': ["C", "c"],
    'T': ["T", "t"],
    'A': ["A", "a"],
    'G': ["G", "g"]
}

# START =========== process excel data =============

def process_sequences(): #feature to process original data
    df = pd.read_csv(r"../Data/CSV_PR_original.csv") 

    # Remove rows with ~ in sequences and rows identical to wildtype sequence
    df = df[~df.NASeq.str.contains('~')]
    df = df[~df.NASeq.str.contains(wild)] 

    # Retrieve same len sequences, drop duplicates and create new df with reset index
    mask = df.NASeq.map(len) == len(wild) 
    #df.drop_duplicates(subset ='NAseq', keep="first", inplace=True)
    df = df[mask].reset_index() 
    return df.NASeq

def process_data(raw = False): 
    seq_data = process_sequences()

    if raw == False: #if raw sequence is false
        seq_data.columns = [seq_col] 
        mask = seq_data[seq_col].map(len) == len(wild) #both aa sequence same len if not remove
        seq_data = seq_data[mask].reset_index()    

    c = [] #empty list
    for seq in seq_data[seq_col]:
        c.append(distance.hamming(wild, seq)) #calculate distance and add to end of empty list c

    c = pd.DataFrame(c) #make list c a df
    
    df = pd.concat([c[0], seq_data], axis=1) #create new df by combining hamming distance and aa sequence
    df.columns = [sub_col, 'index', seq_col] 
    df = df[['index', sub_col, seq_col]] #re order columns
    df.drop_duplicates([seq_col], keep="first", inplace=True) #drop duplicates
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
    for i in range(0, len(sequence), 1):
        codon = sequence[i:i+1]
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

# START =========== grouping ==============

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
    
# END =========== grouping ==============

# START ================= A* ====================

def calc_h(node1, node2): #heuristic function
    n = 0 #number of 0 mutations
    p = 1
    muts = [[x,y] for x,y in zip(node1.sequence, node2.sequence) if x != y]
    for m in muts:
        mp = probs.loc[m[0], m[1]]
        if mp == 0: #mutation = 0
            n = n + 1
        else:
            p = p * mp
    return (1/p) + n

class PriorityQueue:
  def __init__(self):
    self.pq = []

  def print(self):
      for n in self.pq:
          print(n.sequence)

  def add(self, item):
    heappush(self.pq, item)

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
        self.hcost = 0
        self.gcost = self.g()

    def get_neighbours(self):
        for key in list(keys)[self.depth+1:]:
            bin = bins.get_group(key)#.reset_index()
            neighbours = [Node(seq, self.depth+1, self) for seq in bin[seq_col]]
            return neighbours

    def f(self):
        return  self.gcost + self.hcost

    def g(self):
        gcost = 0
        if self.parent != None: #no parent, gcost = 0
            gcost = self.parent.gcost + calc_h(self.parent, self) #got parent, take the parent gcost and calc cost via prob
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
        open.add(start)
        closed = set()
        while open:
            current = open.poll()

            if current.sequence == goal.sequence:
                print("----PATH----")
                path = self.rebuild_path(current)
                for node in reversed(path):
                    print(node.sequence)
                    print()
                break
            neighbours = current.get_neighbours()
            
            n = 0
            for node in neighbours:
                node.hcost = calc_h(node, goal)
                if node not in closed:
                    open.add(node)
                n = n + 1
                if n >= max_n:
                    break
            closed.add(current)

# END ================= A* ====================

data = get_aa_data("data_test.csv") 

probs = get_probs().transpose() #row become column, column become row
bins = data.groupby(sub_col)

print_groups(bins)
keys = bins.groups.keys()
#print(keys)

start = Node(wild, 0, None)
goal = Node(get_seq_from_group(bins, 7, 0), 7, None)

start_time = time.process_time()

# dict_keys([4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,29,30,31,33,34])
depth = 5
index = 0
parent = None
sequence = get_seq_from_group(bins, depth, index)
goal = Node(sequence, depth, parent)

print("-----------")
print("Start: " + str(start.sequence))
print()
print("End: " + str(goal.sequence))

if goal:
    astar = Solver()
    astar.solve(start, goal)

print ("Execution Time: " + str(time.process_time() - start_time))