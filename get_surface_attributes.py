from ase.io import read
import numpy as np
from chem_env_PDH_surf import *
from networkx import *
from sys import argv
from ase.visualize import view
import itertools


def find_d_band(atom,nedos,atoms_dos,E_cutoff):
    print(atom)
    band_index = []
    band_weight = []
    av = atoms_dos[atom].dos['d']["av"]
    for nE in range(nedos):
        if atoms_dos[atom].E[nE] < E_cutoff:
            band_index.append(atoms_dos[atom].E[nE])
            band_weight.append(av[nE])

    center = Moment(band_index, band_weight)
    width = Moment(band_index, band_weight,2)**0.5
    return center, width

def Moment(xdata,ydata,nth=1):
    """
    Internal function to calculate the nth moment  of the
    xdata,ydata distribution, by default the first moment
    (center of energy)
    
    John Kitchin <jkitchin@andrew.cmu.edu>
    """
    xdata=np.array(xdata)
    ydata=np.array(ydata)
    data=(xdata**(nth))*ydata
    normalization=np.trapz(ydata,x=xdata)

    if nth == 0:
        return normalization
    elif normalization > 1.0e-6:
        return np.trapz(data,x=xdata)/normalization
    else:
        return 0



def calculate_dband(path,surf_atoms):
    E_cutoff = 5
    atoms = read(path+'/POSCAR')
    d_file = open(path+'/DOSCAR')
    lines = d_file.readlines()
    chemical_symbols = atoms.get_chemical_symbols()
    natoms = len(chemical_symbols)
    header = lines[:6]  ## These are the lines right before the energies along with occupation starts to print.
    header[5].split()
    e_max = float(header[5].split()[0])-0.1
    e_min = float(header[5].split()[1])+0.1
    nedos = int(header[5].split()[2])
    print('Number of DOS points are: {}'.format(nedos))
    efermi = float(header[5].split()[3])
    print('the fermi energy is:{}'.format(efermi))
    sumdos = lines[6:(nedos+6)] # The total DOS and integrated DOS for the overall system
    split_dos = lines[(nedos+6):] # This is the individual dos for each atoms from here now!
    Ldos = []

    for atom in range(natoms):
        start = atom * (nedos + 1) + 1
        Ldos.append(split_dos[start:start+nedos])

    atoms_dos = []
    
    class AtomDOS(object):
        def __init__(self, symbol, Ldos):
            self.symbol = symbol
            self.E = []
            self.dos = {"s":{"up":[],"dn":[],"av":[]},
                        "p":{"up":[],"dn":[],"av":[]},
                        "d":{"up":[],"dn":[],"av":[]},}

            for nE in range(nedos):
                Ldos_tmp = list(map(float, Ldos[nE].split()))
                self.E.append(Ldos_tmp[0]-efermi)
                
                def spin_decomp(band, up, dn):
                    self.dos[band]["up"].append(Ldos_tmp[up])
                    self.dos[band]["dn"].append(Ldos_tmp[dn])
                    self.dos[band]["av"].append((Ldos_tmp[up] + Ldos_tmp[dn]) / 2)

                spin_decomp("s", 1, 2)
                spin_decomp("p", 3, 4)
                spin_decomp("d", 5, 6)

    #Get the dos states for individual atoms in an array object!

    print(natoms)
    for atom in range(natoms):
        atoms_dos.append(AtomDOS(chemical_symbols[atom], Ldos[atom] ) )
    print('Summing over the given atoms')
    up_s = [0] * nedos
    dn_s = [0] * nedos
    av_s = [0] * nedos
    p_O_up = 0
    p_O_down = 0
    dband_center_avg = []
    dband_width_avg = []
    for b in ['d']:
        for sumi in surf_atoms:
            center,width = find_d_band(int(sumi),nedos,atoms_dos,E_cutoff)
            dband_center_avg.append(center)
            dband_width_avg.append(width)
            for nE in range(nedos):
                up_s[nE] += atoms_dos[sumi].dos[b]["up"][nE]
                dn_s[nE] += atoms_dos[sumi].dos[b]["dn"][nE]
                av_s[nE] += atoms_dos[sumi].dos[b]["av"][nE]
    dband_c_avg = np.mean(dband_center_avg)
    dband_w_avg = np.mean(dband_width_avg)
    sums = surf_atoms
    up_s[:] = [x / len(sums) for x in up_s]
    dn_s[:] = [x / len(sums) for x in dn_s]
    av_s[:] = [x / len(sums) for x in av_s]
    #f = open('d_band_info','a')
    #f.write('{},{}'.format(dband_c_avg,dband_w_avg))
    #f.write('\n')
    #f.close()

    print('The average center and width for the given atoms is:{} and {}'.format(str(round(dband_c_avg,2)),str(round(dband_w_avg,2))))
    return(dband_c_avg,dband_w_avg)

def add_to_list(graph,g_node):
    bond_distance = []
    x = []
    index = []
    for i in graph.neighbors(g_node):
        if 'ads' not in graph.nodes[i]:
            x.append(i)
            bond_distance.append(graph.edges[g_node,str(i)]['dist_edge'])
            index.append(graph.nodes[i]['index'])
    return x,bond_distance,index


def get_ads_list(graph):
    ads_nodes = []
    for n,data in graph.nodes(data='ads'):
        if data == True:
            ads_nodes.append(str(n))
    return ads_nodes


radii_multiplier = 1.1
skin = 0.25
f = open(argv[1])
out = open('surface_info','w')
out.write('Surface_atoms' + '\t' + 'coordination' + '\t' + 'avg_site_bond_dis' + '\t' + 'dband_center' + '\t' + 'dband_width'+ '\n')
for i in f.readlines():
    path = i.replace('\n','')
    atoms = read(path+'/POSCAR')
   # path = i.replace('/' + argv[1].split('/')[-1],'') 
    print(path)
    #view(atoms)
    adsorbate_elements = ['C','H']
    adsorbate_atoms = [atom.index for atom in atoms if atom.symbol in adsorbate_elements]
    nl = NeighborList(natural_cutoffs(atoms, radii_multiplier), self_interaction=False,bothways=True, skin=skin)
    nl.update(atoms)
    full, chem_env = process_atoms(atoms,nl,adsorbate_atoms,radius=2)
    ads_nodes = get_ads_list(chem_env[0])
    nearest_neighbor = []
    bond_distance = []
    dos_atoms = []
    nn_index = []
    for i in ads_nodes:
        nn,bond_distances,index = add_to_list(chem_env[0],i)
        nearest_neighbor.append(nn)
        bond_distance.append(bond_distances)
        nn_index.append(index)
    
    nearest_neighbor=(list(itertools.chain.from_iterable(nearest_neighbor)))
    bond_distance=(list(itertools.chain.from_iterable(bond_distance)))
    nn_index=(list(itertools.chain.from_iterable(nn_index)))
    print(nearest_neighbor)
    print(bond_distance)
    print(nn_index)
    if len(nn_index) != 0:
        dband_c_avg,dband_w_avg=calculate_dband(path,nn_index)
        site_indices = str(nn_index)
        out.write('{}  \t  {}  \t  {}  \t  {}  \t  {}  \n'.format(site_indices,len(nn_index),np.mean(bond_distance),dband_c_avg,dband_w_avg))
    else:
        print('Molecule is physorbed or something else is off, please check path')
        out.write('The molecule is physisorbed or something is off chech this path \n')
'''
next_nearest = []
bond_distance =[]
nn_index = []
for i in nearest_neighbor:
    _,nn,bond_distances = add_to_list(chem_env[0],i)
    next_nearest.append(nn)
    bond_distance.append(bond_distances)
    #nn_index.append(index)
print(next_nearest)
print(bond_distance)
CN_ads = len(nearest_neighbor)
CN_site = np.mean(list(len(i) for i in next_nearest))
print(CN_ads, CN_site)
'''
