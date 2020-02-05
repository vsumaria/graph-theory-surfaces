from ase.io import read
import numpy as np
from chem_env_PDH_surf import *
from networkx import *
from sys import argv
from ase.visualize import view
import itertools


def add_to_list(graph,g_node):
    bond_distance = []
    x = []
    for i in graph.neighbors(g_node):
        if 'ads' not in graph.nodes[i]:
            x.append(i)
            bond_distance.append(graph.edges[g_node,str(i)]['dist_edge'])
            
    return x,bond_distance


def get_ads_list(graph):
    ads_nodes = []
    for n,data in graph.nodes(data='ads'):
        if data == True:
            ads_nodes.append(str(n))
    return ads_nodes


radii_multiplier = 1.1
skin = 0.25

atoms = read(argv[1])
#view(atoms)
adsorbate_elements = ['N','O']
adsorbate_atoms = [atom.index for atom in atoms if atom.symbol in adsorbate_elements]
nl = NeighborList(natural_cutoffs(atoms, radii_multiplier), self_interaction=False,bothways=True, skin=skin)
nl.update(atoms)
full, chem_env = process_atoms(atoms,nl,adsorbate_atoms,radius=2)
ads_nodes = get_ads_list(chem_env[0])
nearest_neighbor = []
bond_distance = []
for i in ads_nodes:
    nn,bond_distances = add_to_list(chem_env[0],i)
    nearest_neighbor.append(nn)
    bond_distance.append(bond_distances)

nearest_neighbor=(list(itertools.chain.from_iterable(nearest_neighbor)))
bond_distance=(list(itertools.chain.from_iterable(bond_distance)))
print(nearest_neighbor)
print(bond_distance)
next_nearest = []
bond_distance =[]
for i in nearest_neighbor:
    nn,bond_distances = add_to_list(chem_env[0],i)
    next_nearest.append(nn)
    bond_distance.append(bond_distances)
print(next_nearest)
print(bond_distance)
CN_ads = len(nearest_neighbor)
CN_site = np.mean(list(len(i) for i in next_nearest))
print(CN_ads, CN_site)
