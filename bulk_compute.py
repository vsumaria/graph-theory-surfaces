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

def make_single_point_job(path_dir,i,atoms):
    template_path = '/depot/jgreeley/data/Pushkar_Sid_Alloy/template'
    sin_pnt_path = '{}/single_point_{}'.format(path_dir,i)
    os.makedirs(sin_pnt_path)
    
    cell_array = [atoms.get_cell()[0][0],atoms.get_cell()[1][1],atoms.get_cell()[2][2]]
    y = atoms.copy()
    for i, cell_len in enumerate(cell_array):
        if cell_len < 6.5 and i ==0:
            y = y.repeat([2,1,1])
        if cell_len < 6.5 and i ==1:
            y = y.repeat([1,2,1])
        if cell_len < 6.5 and i ==2:
            y = y.repeat([1,1,2])
            
    del y[i]
    s = AseAtomsAdaptor.get_structure(read(path_dir+'/'+'CONTCAR'))
    vol_k = int(24**3/s.volume)
    k=Kpoints.automatic_density_by_vol(structure=s,kppvol=vol_k)
    k.write_file(sin_pnt_path+'/'+'KPOINTS')
    y.write('{}/POSCAR'.format(sin_pnt_path))
    shutil.copyfile('{}/POTCAR'.format(path_dir), '{}/POTCAR'.format(sin_pnt_path))
    shutil.copyfile('{}/INCAR'.format(template_path), '{}/INCAR'.format(sin_pnt_path))
    shutil.copyfile('{}/qs_vasp'.format(template_path), '{}/qs_vasp'.format(sin_pnt_path))

for di in os.listdir(argv[1]):
    path_dir = '{}/{}'.format(argv[1],di)
    x = read('{}/CONTCAR'.format(path_dir))
    print(x)

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
            make_single_point_job(path_dir,i,x)
