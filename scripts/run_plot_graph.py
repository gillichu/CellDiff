import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import genfromtxt

mat = np.zeros((7, 7))
filename = "optimal_mat_12h.txt"
with open(filename, "r") as r:
    lines = r.readlines()
    for i, line in enumerate(lines):
        for j, val in enumerate(line.split('\t')):
            mat[i][j] = val
print(mat)
#mat = genfromtxt(filename, delimiter='\t')
#print(mat)

mapping = {0: 'Pluripotent-progenitor',
           1: 'PCGLC',
           2: 'Germlayer-progenitor',
           3: 'NMPs',
           4: 'Endoderm',
           5: 'Neural',
           6: 'Somite'}

G = nx.DiGraph()
# add nodes
for node in mapping.keys():
    G.add_node(node)
# add edges
for i in mat:
    for j in 
G.add_edges_from(mat, directed=True)
nx.draw_networkx(G, with_labels = True, arrows=True, node_color = "c", edge_color = "k", font_size = 8)

plt.axis('off')
plt.draw()
plt.savefig("celltype_summary.pdf")
