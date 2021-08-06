from numpy import *
from scipy import sparse
from scipy import linalg
from itertools import count
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import networkx as nx
import subprocess as subprocess

# PDB file i/o
PDB_ID = sys.argv[1]
# PDB_ID = '1URR'
file_URL = "https://files.rcsb.org/download/%s.pdb" % (PDB_ID)
local_PDB_file = PDB_ID + ".pdb"
local_contacts_file = PDB_ID + "_contacts.tsv"

print(subprocess.run(["wget", file_URL, ">", local_PDB_file], capture_output=True))
print(subprocess.run(["python3", "get_static_contacts.py" --structure $local_PDB_file --itypes all --output local_contacts_file], capture_output=True))python3 get_static_contacts.py --structure $local_PDB_file --itypes all --output local_contacts_file

# ! wget $file_URL > $local_PDB_file
# ! python3 get_static_contacts.py --structure $local_PDB_file --itypes all --output local_contacts_file

# import contacts list
file_path = local_contacts_file
contacts = pd.read_csv(file_path, delimiter = '\t', names = [1, 'i_type', 'node_1', 'node_2'])
contacts_file_name = os.path.basename(file_path)

# generate adjacency matrix
def generate_adj_matrix(contacts):
    vertices = contacts.iloc[:, 2::]
    ea = array(vertices)
    e_atoms = [item for sublist in ea for item in sublist]
    unique_atoms = set(e_atoms)
    numverts = len(unique_atoms)
    adj_matrix = sparse.lil_matrix((numverts, numverts))
    nums = range(numverts)

    atoms_dict = dict(zip(unique_atoms, nums))

    for contact in ea:
        adj_matrix[atoms_dict[contact[0]], atoms_dict[contact[1]]] = 1
        adj_matrix[atoms_dict[contact[1]], atoms_dict[contact[0]]] = 1

    return adj_matrix, atoms_dict

adj_matrix, atoms_dict = generate_adj_matrix(contacts)

# visualize adjacency matrix
def plot_adj_matrix(adj_mat):
    plot_title = "%s Contact Adjacency Matrix" % (PDB_ID)
    plt.spy(adj_mat, markersize = 2)
    plt.ylabel("Atom Index")
    plt.xlabel("Atom Index")
    plt.title(plot_title)
    plt.show()

plot_adj_matrix(adj_matrix)

# generate graph from adjacency matrix
def show_graph_no_labels(adj_matrix, mylabels=None):
    rows, cols = sparse.find(adj_matrix)[0:2]
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    gr.add_nodes_from(nodes)
    nx.draw(gr, node_size=3) # labels=mylabels, with_labels=True)
    plt.savefig(("%s_contact_graph.png") % (PDB_ID))
    plt.show()

    return gr

gr = show_graph_no_labels(adj_matrix,  atoms_dict.items())

# color graph by amino acid index
# groups = set(nx.get_node_attributes(gr,'group').values())
# mapping = dict(zip(sorted(groups),count()))
# nodes = gr.nodes(data=True)
# colors = [mapping[gr.nodes[n]['group']] for n in nodes]
# gr.nodes[0][]
# pos = nx.spring_layout(gr)
# e_col = nx.draw_networkx_edges(gr, pos, alpha=0.2)
# n_col = nx.draw_networkx_nodes(gr, nodelist=nodes, node_color=colors, with_labels=False, node_size=10, cmap=plt.cm.jet)
# nx.draw(gr, node_color=color_map, node_size=20)
# plt.show()

# get graph laplacian


# find conductance
