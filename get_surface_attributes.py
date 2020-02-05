from ase.io import read
import numpy as np
from chem_env_sid import *
from networkx import *
from sys import argv
from ase.visualize import view
import itertools

def add_to_list(graph,g_node):
    x = []
    for i in graph.neighbors(g_node):
        if 'ads' not in graph.nodes[i]:
            x.append(i)

    return x


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
for i in ads_nodes:
    nn = add_to_list(chem_env[0],i)
    nearest_neighbor.append(nn)

nearest_neighbor=(list(itertools.chain.from_iterable(nearest_neighbor)))
print(nearest_neighbor)

next_nearest = []
for i in nearest_neighbor:
    nn = add_to_list(chem_env[0],i)
    next_nearest.append(nn)
print(next_nearest)
#next_nearest_neighbor = add_to_list(chem_env[0],nearest_neighbor)

#print(ads_nodes,nearest_neighbor,next_nearest_neighbor)
