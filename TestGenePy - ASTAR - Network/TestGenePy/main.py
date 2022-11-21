import data_processing as dp
import data_graphing as dgr
import networkx as nx
import matplotlib.pyplot as plt

#n = (10, 50, 100, 1000, 13405)
n = 100
data = dp.get_aa_data("data_aa.csv").head(n)
probs = dp.get_probs()

G = dgr.get_graph("../Data/graphs/graph_" + str(n) + ".txt", data, probs, False, 3)
print(nx.info(G))
print()

components = list(nx.connected_components(G))
subs = (G.subgraph(c).copy() for c in sorted(nx.connected_components(G), key=len, reverse=True))

nv = list(sorted(G.nodes()))

#print all possible paths from one node to another
start = 0
end = 0
for path in nx.all_simple_paths(G, nv[start], nv[end]):
    print(path)

#print full sequence based on the index
seq_id = 0
val = data['aa_seq'].values[seq_id]
print()
print("The sequence in index "+ str(seq_id) + " is " + val)

f = 1
for C in subs:
    pos = nx.spring_layout(C)
    labels_pos = {}
    x_off = 0.025 
    y_off = 0.05 
    for k, v in pos.items():
        labels_pos[k] = (v[0]+x_off, v[1]+y_off)
    plt.figure(f)

    nx.draw(C, pos=pos, node_size=30, with_labels=False)
    nx.draw_networkx_labels(C, pos=labels_pos, font_size=10)
    nx.draw_networkx_edges(C,pos=pos)
    
    plt.show()
    f = f + 1



