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
import matplotlib.pyplot as plt


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

def get_neighbor_pairs(sub_g1,node_c):
    Pd_Pd = 0
    Pd_Sn = 0
    Sn_Sn = 0
    for i in sub_graph.edges():
        if node_c not in i:
            if ((i[0].split(':')[0]) == 'Pd' and (i[1].split(':')[0]) == 'Pd'):
                Pd_Pd += 1
            elif ((i[0].split(':')[0]) == 'Sn' and (i[1].split(':')[0]) == 'Sn'):
                Sn_Sn += 1
            else:
                Pd_Sn +=1
    return Pd_Pd, Pd_Sn, Sn_Sn


def get_env_2(sub_g1,sub_g2):
    Pdn = 0
    Snn = 0
    for n in sub_g2.nodes():
        if n not in list(sub_g1.nodes()):
        #print(n)
            if n.split(':')[0] == 'Pd':
                Pdn+=1
            if n.split(':')[0] == 'Sn':
                Snn+=1
    return Pdn,Snn

bulk_formation = {'Pt':-6.06,'Pd':-5.16,'Sn':-3.83}
print(bulk_formation['Pt'])
mix = 0.25
out = open('vacancy_similarity_info','w')
out.write('Path' + '\t' + 'Unique_atom' + '\t' + 'Unique_atom_type' + '\t' + 'similarity_radius_1' + '\t' + 'similarity_radius_2'+ '\t'+ '\t' + 'similarity_radius_total' +  '\t' +  'Pd_Pd_pairs'  +  '\t' +  'Pd_Sn_pairs' +  '\t' +  'Sn_Sn_pairs'+ '\t'  + 'vacancy_formation' + '\t'  + 'total_energy'+ '\n')

Pd_conc = []
Pd_conc_array = []

for di in os.listdir(argv[1]):
    path_dir = '{}/{}'.format(argv[1],di)
    try:
        x = read('{}/OUTCAR'.format(path_dir))
        print(path_dir)
    except:
        print('cannot read here {}'.format(path_dir))
        continue


    nl = NeighborList(natural_cutoffs(x, 1.1), self_interaction=False, bothways = True,skin=0.25)
    nl.update(x)

    full,chem = process_atoms(x, nl, adsorbate_atoms=[], radius=2, grid_n=[1,1,1])
    unique_sub_graphs = []
    unique_index = []
    for i in range(0,len(x)):
        node_c = '{}:{}[0,0,0]'.format(x[i].symbol,i)
        sub_graph = nx.ego_graph(full,node_c,radius =2)
        for j,uni in enumerate(unique_sub_graphs):
            if compare_chem_envs([sub_graph],[uni]):     ## the bond attribute of the graph edges are being compared here and the symb\ols in the edges matter.
                print('atom {} chemical environment is similar to atom {}'.format(i,j))
                break
        else:
            print('this atom is unique: {}, and has a symbol: {}'.format(i,x[i].symbol))
            unique_sub_graphs.append(sub_graph)
            unique_index.append(i)
            Pdn,Snn, sub_graph1 = get_env(full,node_c)
            Pd_1 = Pdn/(Pdn+Snn)
            Sn_1 = Snn/(Pdn+Snn)
            print('First shell contribution:')
            print(Pd_1,Sn_1)

            ##### calculate second shell contribution here #############

            sub_graph2 = nx.ego_graph(full,node_c,radius =2)
            sub_graph2.remove_node(node_c) ## remove self_node
            
            Pdn2,Snn2 = get_env_2(sub_graph1,sub_graph2)
            Pd_2 = Pdn2/(Pdn2+Snn2)
            Sn_2 = Snn2/(Pdn2+Snn2)
            
            print('Directly from graph analysis, second coordination shell is')
            print(Pd_2,Sn_2)

            Pdo = (mix*np.mean(Pd_2)+Pd_1)/(1+mix)
            Sno = (mix*np.mean(Sn_2)+Sn_1)/(1+mix)
            print('First + second contribution is')
            print(Pdo,Sno)
            Pd_Pd, Pd_Sn, Sn_Sn = get_neighbor_pairs(sub_graph1,node_c)
            t_bi = Pd_Pd + Pd_Sn + Sn_Sn
            Pd_Pd_r = Pd_Pd/t_bi
            Sn_Sn_r = Sn_Sn/t_bi
            Pd_Sn_r = Pd_Sn/t_bi
            path_vac = path_dir + '/' + 'single_point_{}'.format(i)
            outcar_vac = read(path_vac+'/'+'OUTCAR')
            vac_energy = outcar_vac.get_potential_energy()
            rep = int(outcar_vac.get_number_of_atoms()/x.get_number_of_atoms())+1
            tot_ener = x.get_potential_energy() * rep
            vac_formation = vac_energy + bulk_formation[x[i].symbol] - tot_ener
            
            ######## Make Histograms here #############
            if Pd_1 in Pd_conc:
                index = Pd_conc.index(Pd_1)
                print(index)
                Pd_conc_array.append([])
                Pd_conc_array[index].append(vac_formation)
            else:
                Pd_conc.append(Pd_1)
                index = Pd_conc.index(Pd_1)
                print(index)
                Pd_conc_array.append([])
                Pd_conc_array[index].append(vac_formation)

            print(vac_formation)
            out.write(path_dir + '\t' + str(i) + '\t' + x[i].symbol + '\t' + '{},{}'.format(Pd_1,Sn_1) + '\t' + '{},{}'.format(Pd_2,Sn_2) +  '\t' + '{},{}'.format(Pdo,Sno)  +  '\t' +  '{}'.format(Pd_Pd_r)  +  '\t' +  '{}'.format(Pd_Sn_r) +  '\t' +  '{}'.format(Sn_Sn_r) +  '\t' +  '{}'.format(vac_formation)  +  '\t' +  '{}'.format(tot_ener/rep) + '\n')

print('Plotting histograms now')
for i,n  in enumerate(Pd_conc):
    plt.figure()
    binwidth = 0.05
    data = Pd_conc_array[i]
    N = (max(data)-min(data))/binwidth + 1
    binlist = np.linspace(min(data),max(data),N)
    plt.hist(data, bins=binlist)
    #plt.hist(Pd_conc_array[i])
    plt.savefig('Pd_{}_his.png'.format(round(n,2)),dpi=400)

#### Below is an alternative method to get second shell contribution.

'''
Pds = []
Sns = []

for node in sub_graph1.nodes():
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
'''
