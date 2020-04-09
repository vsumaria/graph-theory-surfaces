#! /usr/bin/env python

from ase.io import read,write
import numpy as np
from chem_env_PDH_surf import *
from networkx import *
from sys import argv
from ase.visualize import view
import itertools
import os
import shutil
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.ase import AseAtomsAdaptor
from os import system



def get_env(full,node_c):

    
    sub_graph = nx.ego_graph(full,node_c,radius =1)
    sub_graph.remove_node(node_c) ## remove self_node

    Pdn = 0
    Snn = 0
    for n in sub_graph.nodes():
        #print(n)
        if n.split(':')[0] == 'Pd':
            Pdn+=1
        if n.split(':')[0] == 'Sn':
            Snn+=1
    return Pdn,Snn, sub_graph

x = read('test/POSCAR_PdSn_3_5')
i = 1  ### specify the atom around which field of composition has to be estimated

node_c = '{}:{}[0,0,0]'.format(x[i].symbol,i)

nl = NeighborList(natural_cutoffs(x, 1.1), self_interaction=False, bothways = True,skin=0.25)
nl.update(x)

full,chem = process_atoms(x, nl, adsorbate_atoms=[], radius=2, grid_n=[1,1,1])
Pdn,Snn, sub_graph = get_env(full,node_c)
Pd_1 = Pdn/(Pdn+Snn)
Sn_1 = Snn/(Pdn+Snn)
print('First shell contribution:')
print(Pd_1,Sn_1)

##### calculate second shell contribution here #############

Pds = []
Sns = []

for node in sub_graph.nodes():
    print('calculating env. for node: {}'.format(node))
    Pdn,Snn, sub_graph = get_env(full,node)
    Pd_1_s = Pdn/(Pdn+Snn)
    Sn_1_s = Snn/(Pdn+Snn)
    print(Pd_1_s,Sn_1_s)
    Pds.append(Pd_1_s)
    Sns.append(Sn_1_s)

print('mean second shell contributions are:')
print(np.mean(Pds),np.mean(Sns))

mix = 0.25

print('First + second contribution is')

Pdo = (mix*np.mean(Pds)+Pd_1)/(1+mix)
Sno = (mix*np.mean(Sns)+Sn_1)/(1+mix)
print(Pdo,Sno)