# Identifies atomic contacts within a pdb structure using GetContacts
# framework (courtesy of Dror Lab) and visualizes them as an unweighted,
# undirected graph of these contacts.

import numpy as np
import scipy
from scipy import sparse
import itertools
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import networkx as nx
import subprocess
import Bio
from Bio.PDB.PDBParser import PDBParser

# name PDB file variables
# PDB_ID = '1URR'
PDB_ID = sys.argv[1]
store_directory = './structures/%s' % PDB_ID
if not os.path.exists(store_directory):
    os.makedirs(store_directory)
if not os.path.exists(store_directory):
    os.makedirs(store_directory)
file_URL = "https://files.rcsb.org/download/%s.pdb" % (PDB_ID)
local_PDB_file = "structures/%s/%s.pdb" % (PDB_ID, PDB_ID)
local_contacts_file = "structures/%s/%s_contacts.tsv" % (PDB_ID, PDB_ID)

# download pdb file and generate contacts tsv
bashCommand = "wget -O %s %s" % (local_PDB_file, file_URL)
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

bashCommand = "get_static_contacts.py --structure %s --itypes all --output %s" % (local_PDB_file, local_contacts_file)
process = subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)

# import contacts list
contacts = pd.read_csv(local_contacts_file, delimiter = '\t', skiprows=[0, 1], names = ['frame', 'i_type', 'node_1', 'node_2'])

# generate adjacency matrix
def generate_adj_matrix(contacts):
    vertices = contacts.iloc[:, 2::]
    ea = np.array(vertices)
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

# generate atom dictionaries
adj_matrix, atoms_dict = generate_adj_matrix(contacts)
atom_index_dict = {v: dict(zip(["frame", "res_name", "res_idx", "atom_ID"], k.split(":"))) for k, v in atoms_dict.items()}

# visualize adjacency matrix
def plot_adj_matrix(adj_mat):
    plot_title = "%s Contact Adjacency Matrix" % (PDB_ID)
    plt.spy(adj_mat, markersize = 0.5)
    plt.ylabel("Atom Index")
    plt.xlabel("Atom Index")
    plt.title(plot_title)
    plt.show()

plot_adj_matrix(adj_matrix)

# generate graph from adjacency matrix
def show_graph_no_labels(adj_matrix, mylabels=None):
    rows, cols = scipy.sparse.find(adj_matrix)[0:2]
    edges = zip(rows.tolist(), cols.tolist())
    plot_title = "%s Contact Network" % (PDB_ID)
    graph = nx.Graph()
    graph.add_edges_from(edges)
    graph.add_nodes_from(atom_index_dict.items())
    pos = nx.spring_layout(graph)
    options = {"node_size": 0.1, "alpha": 0.9}
    nx.draw_networkx_edges(graph, pos, alpha=0.3)
    nx.draw_networkx_nodes(graph, pos, **options)
    plt.title(plot_title)
    plt.savefig(("%s/%s_contact_graph.png") % (PDB_ID, PDB_ID), dpi=300)
    plt.show()

    return graph

graph = show_graph_no_labels(adj_matrix,  atoms_dict.items())

# color graph by residue number
def show_graph_aa_color(graph):
    aa_idxs = set(nx.get_node_attributes(graph,'res_idx').values())
    mapping = dict(zip(sorted(aa_idxs), itertools.count()))
    nodes = graph.nodes()
    colors = [mapping[nodes[n]['res_idx']] for n in nodes]

    pos = nx.spring_layout(graph)
    ec = nx.draw_networkx_edges(graph, pos, alpha=0.2)
    nc = nx.draw_networkx_nodes(graph, pos, nodelist=nodes, node_color=colors, node_size=5, cmap=plt.get_cmap('viridis'))

    plot_title = "%s Colored by Residue Number" % (PDB_ID)
    plt.colorbar(nc)
    plt.axis('off')
    plt.title(plot_title)
    plt.savefig(("%s/%s_contact_graph_aa_idx.png") % (PDB_ID, PDB_ID), dpi=300)
    plt.show()

show_graph_aa_color(graph)

# color graph by residue depth
# parser = PDBParser()
# structure = parser.get_structure(PDB_ID, local_PDB_file)
# model = structure[0]
# rd = Bio.PDB.ResidueDepth(model)
# print(rd['A',(' ', 152, ' ')])

#
# # get eigenvalues of graph laplacian
#
# def get_laplacian(A):
#     n,m = A.shape
#     diags = A.sum(axis=1)
#     D = scipy.sparse.spdiags(diags.flatten(), [0], m, n, format='csr')
#     lap = D - A
#
#     return lap
#
# lap = get_laplacian(adj_matrix)
#
# vals, vecs = scipy.sparse.linalg.eigs(lap)
#
# sorted(vecs[-2]
#
#
# vals
# vecs[:, 0]
#
# # other algorithms to identify communities in network
# G = graph
# all_connected_subgraphs = []
#
# graphs = sorted(nx.connected_components(graph), key=lambda x: -len(x))
# largest_cc = graph.subgraph(graphs[0])
#
# nx.draw(largest_cc, node_size=1)
#
#
# comms = nx.algorithms.community.asyn_fluidc(largest_cc, 5)
# list(comms)[1]
#
# G = nx.path_graph(4)
# G.add_edge(5,6)
# graphs = sorted(nx.connected_components(G))
# # calculate phi-psi angles
# # for model in Bio.PDB.PDBParser().get_structure("1URR", "1URR.pdb.1"):
# #     for chain in model :
# #         poly = Bio.PDB.Polypeptide.Polypeptide(chain)
# #         print "Model %s Chain %s" % (str(model.id), str(chain.id)),
# #         print poly.get_phi_psi_list()
#
# # find conductance
