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

def K_POINT_make(path,repeat_vector,sin_pnt_path):
    kp = []
    #print(repeat_vector)
    for i,rep in enumerate(repeat_vector):
        f = open(path+'/KPOINTS')
        kp.append(int(float(f.readlines()[3].split(' ')[i])/rep))
    #print(kp)
    new_f = open(sin_pnt_path+'/KPOINTS','w')
    f = open(path+'/KPOINTS')
    for i, kps in enumerate(f.readlines()):
        if i != 3:
            new_f.write(kps)
        else:
            for k in kp:
                new_f.write('{} '.format(str(k)))
            new_f.write('\n')
    f.close()
    new_f.close()

def make_single_point_job(path_dir,i,atoms):
    template_path = '/depot/jgreeley/data/Pushkar_Sid_Alloy/template_single'
    sin_pnt_path = '{}/single_point_{}'.format(path_dir,i)
    os.makedirs(sin_pnt_path)
    
    cell_array = [atoms.get_cell()[0][0],atoms.get_cell()[1][1],atoms.get_cell()[2][2]]
    y = atoms.copy()
    tot_atoms = len(atoms)
    cell_cutoff = 7.5
    repeat_array = [1,1,1]
    for j, cell_len in enumerate(cell_array):
        if cell_len < cell_cutoff and j ==0:
            mult = int(cell_cutoff/cell_len)+1
            y = y.repeat([mult,1,1])
            repeat_array[j] = 2 
            tot_atoms = tot_atoms*mult
        if cell_len < cell_cutoff and j ==1:
            mult = int(cell_cutoff/cell_len)+1
            y = y.repeat([1,mult,1])
            tot_atoms =tot_atoms*mult
            repeat_array[j] = 2
        if cell_len < cell_cutoff and j ==2:
            mult = int(cell_cutoff/cell_len)+1
            y = y.repeat([1,1,mult])
            tot_atoms =tot_atoms*mult
            repeat_array[j] = 2
    print(i,y[i].symbol)
    del y[i]
    K_POINT_make(path_dir,repeat_array,sin_pnt_path)
    nodes = int(tot_atoms/20)
    y.write('{}/POSCAR'.format(sin_pnt_path))
    shutil.copyfile('{}/POTCAR'.format(path_dir), '{}/POTCAR'.format(sin_pnt_path))
    shutil.copyfile('{}/INCAR'.format(template_path), '{}/INCAR'.format(sin_pnt_path))
    shutil.copyfile('{}/qs_vasp'.format(template_path), '{}/qs_vasp'.format(sin_pnt_path))
    #s = AseAtomsAdaptor.get_structure(read(sin_pnt_path+'/'+'POSCAR'))
    #vol_k = int(24**3/s.volume)
    #k=Kpoints.automatic_density_by_vol(structure=s,kppvol=vol_k)
    #k.write_file(sin_pnt_path+'/'+'KPOINTS')
    system("cd "+sin_pnt_path+" ; change_nodes.sh {}".format(nodes))
    system("cd "+sin_pnt_path+" ; sortatoms.py POSCAR")
    system("cd "+sin_pnt_path+" ; NIFEmakePOT")
    system("cd "+sin_pnt_path+" ; sbatch qs_vasp")

for di in os.listdir(argv[1]):
    path_dir = '{}/{}'.format(argv[1],di)
    try:
        x = read('{}/CONTCAR'.format(path_dir))
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
        node_c = '{}:{}[0,0,0]'.format(x[i].symbol,i)    #### choose the center node for each atom in the unit cell
        sub_graph = nx.ego_graph(full,node_c,radius =2)  #### make a subgraph out of that node 
        for j,uni in enumerate(unique_sub_graphs):
            if compare_chem_envs([sub_graph],[uni]):     ## the bond attribute of the graph edges are being compared here and the symbols in the edges matter.
                print('atom {} chemical environment is similar to atom {}'.format(i,j))
                break
        else:
            print('this atom is uniue: {}'.format(i))
            unique_sub_graphs.append(sub_graph)
            unique_index.append(i)
            make_single_point_job(path_dir,i,x)
