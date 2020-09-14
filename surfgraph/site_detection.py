from numpy.linalg import norm
from ase.neighborlist import NeighborList, natural_cutoffs
from itertools import combinations
from ase.constraints import constrained_indices

from surfgraph.chemical_environment import process_atoms
from surfgraph.chemical_environment import process_site
from surfgraph.chemical_environment import unique_chem_envs
from surfgraph.chemical_environment import bond_match
from surfgraph.helpers import draw_atomic_graphs
from surfgraph.helpers import normalize
from surfgraph.helpers import offset_position

import networkx.algorithms.isomorphism as iso
import numpy as np
import networkx as nx


def generate_normals_original(atoms, surface_normal=0.5, normalize_final=True, adsorbate_atoms=[]):
    normals = np.zeros(shape=(len(atoms), 3), dtype=float)

    atoms = atoms.copy()

    del atoms[adsorbate_atoms]

    cutoffs = natural_cutoffs(atoms)

    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    cell = atoms.get_cell()

    for index, atom in enumerate(atoms):
        normal = np.array([0, 0, 0], dtype=float)
        for neighbor, offset in zip(*nl.get_neighbors(index)):
            direction = atom.position - offset_position(atoms, neighbor, offset)
            normal += direction
        if norm(normal) > surface_normal:
            normals[index,:] = normalize(normal) if normalize_final else normal

    surface_mask = [index for index in range(len(atoms)) if norm(normals[index]) > 1e-5]

    return normals, surface_mask

def generate_site_type(atoms, surface_mask, normals, coordination, unallowed_elements=[]):
    cutoffs = natural_cutoffs(atoms)

    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    surface_mask = [index for index in surface_mask if atoms[index].symbol not in unallowed_elements]

    possible = list(combinations(set(surface_mask), coordination))
    valid = []
    sites = []

    for cycle in possible:
       for start, end in combinations(cycle, 2):
           if end not in nl.get_neighbors(start)[0]:
               break
       else: # All were valid
            valid.append(list(cycle))

    for cycle in valid:
        tracked = np.array(atoms[cycle[0]].position, dtype=float)
        known = np.zeros(shape=(coordination, 3), dtype=float)
        known[0] = tracked
        for index, (start, end) in enumerate(zip(cycle[:-1], cycle[1:])):
            for neighbor, offset in zip(*nl.get_neighbors(start)):
                if neighbor == end:
                    tracked += offset_position(atoms, neighbor, offset) - atoms[start].position
                    known[index + 1] = tracked

        average = np.average(known, axis=0)

        normal = np.zeros(3)
        for index in cycle:
            neighbors = len(nl.get_neighbors(index)[0])
            normal += normals[index] * (1/neighbors)
        normal = normalize(normal)
        if coordination ==2:
            average[2] = average[2] - 0.5
        if coordination == 3:
            average[2] = average[2] -0.7
            #print(average)
            #print(average[2])
        site_ads =Site(cycle=cycle, position=average, normal=normal)
        sites.append(site_ads)
        
    return sites

def generate_site_graphs(atoms, full_graph, nl, sites, adsorbate_atoms=[], radius=3):
    cutoffs = natural_cutoffs(atoms)

    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    site_envs = [None] * len(sites)
    for index, site in enumerate(sites):
        new_site = process_site(atoms, full_graph, nl, site.cycle, radius=radius)
        site_envs[index] = [new_site]
        site.graph = site_envs[index]

    unique_envs, unique_sites = unique_chem_envs(site_envs, sites)

    return unique_sites

class Site(object):
    def __init__(self, cycle, position, normal, graph=None):
        self.cycle = cycle
        self.position = position
        self.normal = normal
        self.graph = graph

    def __eq__(self, other):
        return nx.is_isomorphic(self.graph, other.graph, edge_match=bond_match)

    def __repr__(self):
        return "Cycle:{}, Position:[{}, {}, {}], Normal:[{}, {}, {}], Graph:{}".format(
                       self.cycle, 
                       self.position[0], self.position[1], self.position[2],
                       self.normal[0], self.normal[1], self.normal[2],
                       self.graph)

    def adsorb(self, atoms, adsorbate, adsorbate_atoms, height=2):
        ads_copy = adsorbate.copy()
        ads_copy.rotate([0, 0, 1], self.normal, center=[0,0,0])
        ads_copy.translate(self.position + (self.normal*height))
        atoms.extend(ads_copy)

        index_to_check = range(len(atoms)-len(ads_copy), len(atoms))

        dist = float("inf")

        if len(adsorbate_atoms) != 0:
            for index in index_to_check:
                dists = atoms.get_distances(index, adsorbate_atoms, mic=True)
                dist = min(dist, dists.min())

        return dist
