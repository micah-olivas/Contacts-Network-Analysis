from numpy import *
from scipy import sparse
from scipy import linalg
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import networkx as nx

# PDB file i/o
PDB_ID = "1URR"
file_URL = "https://files.rcsb.org/download/%s.pdb" % (PDB_ID)
local_PDB_file = PDB_ID + ".pdb"
local_contacts_file = PDB_ID + "_contacts.tsv"

! wget $file_URL > $local_PDB_file
! python3 get_static_contacts.py --structure $local_PDB_file --itypes all --output local_contacts_file

# import contacts list
file_path = local_contacts_file
contacts = pd.read_csv(file_path, delimiter = '\t', names = [1, 'i_type', 'node_1', 'node_2'])
contacts_file_name = os.path.basename(file_path)
structure = contacts_file_name.split('_')[0]
print("successfully imported %s" % (structure))

# generate adjacency matrix
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

# visualize adjacency matrix
def plot_adj_matrix(adj_mat):
    plt.spy(adj_mat, markersize = 2)
    plt.ylabel("Atom Index")
    plt.xlabel("Atom Index")
    plt.show()

plot_adj_matrix(adj_matrix)

# generate graph from adjacency matrix
def show_graph_with_labels(adj_matrix, mylabels=None):
    rows, cols = sparse.find(adj_matrix)[0:2]
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=2) # labels=mylabels, with_labels=True)
    plt.savefig(("%s_contact_graph.png") % (structure))
    plt.show()

rows, cols = sparse.find(adj_matrix)[0:2]
edges = zip(rows.tolist(), cols.tolist())
gr = nx.Graph()
gr.add_edges_from(edges)
nx.draw(gr, with_labels=False) # with_labels=True)
plt.show()

show_graph_with_labels(adj_matrix,  atoms_dict.items())

# get graph laplacian

# find conductance
