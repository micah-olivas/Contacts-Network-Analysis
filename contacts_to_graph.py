# Identifies atomic contacts within a pdb structure using GetContacts
# framework (courtesy of Dror Lab) and visualizes them as an unweighted,
# undirected graph of these contacts.

import numpy as np
import scipy
from scipy import sparse
from numpy import linalg
import itertools
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import networkx as nx
import subprocess
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file


# name PDB file variables
PDB_ID = sys.argv[1]
if ".csv" in PDB_ID:
    csv = pd.read_csv(PDB_ID, names=["PDB ID"])

store_directory = 'structures/%s' % PDB_ID
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
contacts = pd.read_csv(local_contacts_file, sep = '\t', skiprows=[0, 1], names = ['frame', 'i_type', 'node_1', 'node_2'])
contacts = pd.read_csv('/Users/micaholivas/Desktop/Contacts-Network-Analysis/%s' % local_contacts_file, sep = '\t', skiprows=[0, 1], names = ['frame', 'i_type', 'node_1', 'node_2'])


# generate atom, residue, and dssp dictionaries
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
adj_matrix, atoms_dict = generate_adj_matrix(contacts)

dssp_tuple = dssp_dict_from_pdb_file(local_PDB_file)
dssp_dict = dssp_tuple[0]

atom_index_dict = {v: dict(zip(["frame", "res_name", "res_idx", "atom_ID"], k.split(":"))) for k, v in atoms_dict.items()}
for i in atom_index_dict.values():
    i['res_idx'] = int(i['res_idx'])

# add dssp output to dictionary of atoms
for gr_k, gr_v in atom_index_dict.items():
    for dssp_k, dssp_v in dssp_dict.items():
        if gr_v['res_idx'] == dssp_k[1][1]:
            gr_v['asa'] = dssp_v[2]


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
    plt.savefig(("structures/%s/%s_contact_graph.png") % (PDB_ID, PDB_ID), dpi=300)
    plt.show()

    return graph, pos
graph, pos = show_graph_no_labels(adj_matrix,  atoms_dict.items())

# color graph by residue number
def show_graph_res_order_color(graph, pos):
    aa_idxs = set(nx.get_node_attributes(graph,'res_idx').values())
    mapping = dict(zip(sorted(aa_idxs), itertools.count()))
    nodes = graph.nodes()
    colors = [mapping[nodes[n]['res_idx']] for n in nodes]

    ec = nx.draw_networkx_edges(graph, pos, alpha=0.2)
    nc = nx.draw_networkx_nodes(graph, pos, nodelist=nodes, node_color=colors, node_size=5, cmap=plt.get_cmap('viridis'))

    plot_title = "%s Colored by Residue Number" % (PDB_ID)
    plt.colorbar(nc)
    plt.axis('off')
    plt.title(plot_title)
    plt.savefig(("structures/%s/%s_contact_graph_aa_idx.png") % (PDB_ID, PDB_ID), dpi=300)
    plt.show()
show_graph_res_order_color(graph, pos)

# color graph by residue depth
def show_graph_asa_color(graph, pos):
    aa_idxs = set(nx.get_node_attributes(graph,'asa').values())
    mapping = dict(zip(sorted(aa_idxs), itertools.count()))
    nodes = graph.nodes()
    colors = [mapping[nodes[n]['asa']] for n in nodes]

    ec = nx.draw_networkx_edges(graph, pos, alpha=0.2)
    nc = nx.draw_networkx_nodes(graph, pos, nodelist=nodes, node_color=colors, node_size=5, cmap=plt.get_cmap('plasma'))

    plot_title = "%s Colored by Solvent Accessibility" % (PDB_ID)
    plt.colorbar(nc)
    plt.axis('off')
    plt.title(plot_title)
    plt.savefig(("structures/%s/%s_contact_graph_asa_idx.png") % (PDB_ID, PDB_ID), dpi=300)
    plt.show()
show_graph_asa_color(graph, pos)

# calculate phi-psi angles

# find conductance

# # Visualize structures
# import __main__
# __main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
#
# import sys, time, os
# import pymol
# # Load Structures
# pymol.finish_launching()
#
# # Read User Input
# spath = os.path.abspath('1URR.pdb')
# # spath = "/opt/anaconda3/PyMOL.app"
# sname = spath.split('/')[-1].split('.')[0]
#
# pymol.cmd.load(spath, sname)
# pymol.cmd.disable("all")
# pymol.cmd.enable(sname)
# pymol.cmd.png("my_image.png")
#
# pymol.cmd.quit()
#
#
# # K-means cluster graph and visualize
#
# # get eigenvalues of graph laplacian
# def get_laplacian(A):
#     n,m = A.shape
#     diags = A.sum(axis=1)
#     D = scipy.sparse.spdiags(diags.flatten(), [0], m, n, format='csr')
#     lap = D - A
#
#     return lap
#
#     lap = get_laplacian(adj_matrix)
#
#     vals, vecs = scipy.sparse.linalg.eigs(lap)
#
#     len(vecs)
#     len(vals)
#
#
#     my_least = sorted(vecs[4])
#     plt.plot(my_least)
#
#
#     vals
#     vecs[:, 0]
#     vecs
#
#     plt.plot(my_least)
#
#     # other algorithms to identify communities in network
#     G = graph
#     all_connected_subgraphs = []
#
#     graphs = sorted(nx.connected_components(graph), key=lambda x: -len(x))
#     largest_cc = graph.subgraph(graphs[0])
#
#     nx.draw(largest_cc, node_size=1)
#
#
#     comms = nx.algorithms.community.asyn_fluidc(largest_cc, 5)
#     list(comms)[1]
#
#     G = nx.path_graph(4)
#     G.add_edge(5,6)
#     graphs = sorted(nx.connected_components(G))
