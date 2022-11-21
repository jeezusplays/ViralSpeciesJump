from collections import namedtuple
import networkx as nx
import data_processing as dp
from os import path

Node = namedtuple("Node", ["idx", "parent_idx", "n", "neighbours"])

def get_graph(filename, data, probs, force=False, n_muts=1):
    if path.exists(filename) and force==False:
        print("Found Graph")
        G = nx.read_gpickle(filename)
    else:
        print("Building Graph")
        nodes = create_nodes(data, probs, n_muts)
        G = create_graph(nodes)
        nx.write_gpickle(G, filename)
    return G

def create_graph(nodes):
    G = nx.Graph()
    for node in nodes:
        G.add_node(node.idx)
        for idx in node.neighbours:
            G.add_edge(node.idx, idx)
    return G

def create_nodes(data, probs, n_muts):
    prev_x = 0
    children = []
    gen = ((ix, iy, data.iloc[ix, 2], data.iloc[iy, 2]) for ix in range(len(data)) for iy in range(len(data)))
    for ix, iy, x, y in gen:
        if ix != prev_x:
            node = Node(idx=prev_x, parent_idx=None, n=len(children), neighbours=children)
            children = []
            prev_x = ix
            yield node
        snps = dp.get_SNPs(x, y)
        if 1 <= len(snps) <= n_muts:
            p = calc_prob(x, y, probs)
            if p>0 and ix != iy:
                children.append(iy)
    #For handling the last node
    last = len(data)-1
    children = []
    for iy in range(last):
        snps = dp.get_SNPs(data.iloc[last, 2], data.iloc[iy, 2])
        if 1 <= len(snps) <= n_muts:
            p = calc_prob(data.iloc[last, 2], data.iloc[iy, 2], probs)
            if p>0 and last != iy:
                children.append(iy)
    node = Node(idx=last, parent_idx=None, n=len(children), neighbours=children)
    yield node

def calc_prob(s1, s2, probs):
    n = 0
    p = 1
    muts = [[x,y] for x,y in zip(s1, s2) if x != y]
    for m in muts:
        mp = probs.loc[m[0], m[1]]
        if mp == 0:
            return 0
        else:
            p = p * mp
    return p 