from ase.io import read
import numpy as np
from chem_env_PDH_surf import *
from networkx import *
from sys import argv
from ase.visualize import view
import itertools

x = read('/depot/jgreeley/data/Pushkar_Sid_Alloy/Sn_alloys/Sn92Pd8/1-12/CONTCAR')

nl = NeighborList(natural_cutoffs(x, 1.1), self_interaction=False, bothways = True,skin=0.25)
nl.update(x)

full,chem = process_atoms(x, nl, adsorbate_atoms=[], radius=2, grid_n=[1,1,1])

unique_sub_graphs = []
unique_index = []
for i in range(0,len(x)):
    node_c = '{}:{}[0,0,0]'.format(x[i].symbol,i)    #### choose the center node for each atom in the unit cell
    sub_graph = nx.ego_graph(full,node_c,radius =2)  #### make a subgraph out of that node 
    for j,uni in enumerate(unique_sub_graphs):
        if compare_chem_envs([sub_graph],[uni]):
            print('atom {} chemical environment is similar to atom {}'.format(i,j))
            break
    else:
        print('this atom is uniue: {}'.format(i))
        unique_sub_graphs.append(sub_graph)
        unique_index.append(i)
