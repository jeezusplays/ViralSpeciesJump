from os import path
import scipy.spatial as distance
import pandas as pd

sub_col = 'num_snp'
seq_col = 'aa_seq'

wild = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"

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

def process_sequences(): #feature to process original data
    df = pd.read_csv(r"../Data/CSV_PR_original.csv") 
    # Remove rows with ~ in sequences and rows identical to wildtype sequence
    df = df[~df.NASeq.str.contains('~')]
    df = df[~df.NASeq.str.contains(wild)] 
    # Retrieve same len sequences, drop duplicates and create new df with reset index
    mask = df.NASeq.map(len) == len(wild) 
    df = df[mask].reset_index() 
    return df.NASeq

def process_data(raw = False):
    seq_data = process_sequences()

    if raw == False: #if raw sequence is false
        wild_aa = ''.join(aa_gen(wild)) #create wildtype amino acid
        seq_data = pd.DataFrame(aa_gen(seq_data)) #create amino sequences from data sequences
        seq_data.columns = [seq_col] 
        mask = seq_data[seq_col].map(len) == len(wild_aa) #both aa sequence same len if not remove
        seq_data = seq_data[mask].reset_index()    

    c = [] #empty list
    for seq in seq_data[seq_col]:
        c.append(distance.levenshtein(wild_aa, seq)) #calculate distance and add to end of empty list c
    c = pd.DataFrame(c) #make list c a df
    
    df = pd.concat([c[0], seq_data], axis=1) #create new df by combining hamming distance and aa sequence
    df.columns = [sub_col, 'index', seq_col] 
    df = df[['index', sub_col, seq_col]] #re order columns
    df.drop_duplicates([seq_col], keep="first", inplace=True) #drop duplicates
    return df
  
def get_raw_data(filename):
    if path.exists("../Data/" + filename):
        print("Found raw data file")
        data = pd.read_csv(r"../Data/" + filename)
    else:
        data = process_data(True)
        data.to_csv(r"../Data/" + filename, index=False)
    return data

def get_aa_data(filename):
    if path.exists("../Data/" + filename):
        print("Found data file")
        data = pd.read_csv(r"../Data/" + filename)
    else:
        data = process_data(False)
        data.to_csv(r"../Data/" + filename, index=False)
    return data

def get_probs():
    if path.exists("../Data/probs.csv"):
        print("Found probs file")
        probs = pd.read_csv(r"../Data/probs.csv").set_index('AA')
    else:
        probs = process_SNP_probabilities()
        probs.to_csv(r"../Data/probs.csv")
    return probs

def aa_multi_sequence_gen(seq_data):
    for sequence in seq_data:
        if len(sequence) != len(wild):
            print("aa_multi_sequence_gen() error: len mismatch: " + str(len(sequence)))
        aa_sequence = ''.join(aa_gen(sequence))
        yield aa_sequence

def aa_gen(sequence):
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        aa = codon_aa(codon.upper())
        if aa == None:
            print("aa_gen(): codon error " + codon)
            continue
        else:
            yield aa

def codon_aa(codon):
    for key, value in codons.items():
        if(codon in value):
            return key
