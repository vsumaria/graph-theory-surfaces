from ase.io import read
import numpy as np
from chem_env_PDH_surf import *
from networkx import *
from sys import argv
from ase.visualize import view
import itertools

x = read('/depot/jgreeley/data/Pushkar_Sid_Alloy/Sn_alloys/Pd100/Sn/CONTCAR')

nl = NeighborList(natural_cutoffs(x, 1.1), self_interaction=False, bothways = True,skin=0.25)
nl.update(x)

full,chem = process_atoms(x, nl, adsorbate_atoms=[], radius=1, grid_n=[1,1,1])
