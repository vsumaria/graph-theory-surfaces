from ase.data.colors import jmol_colors
from ase.neighborlist import NeighborList, natural_cutoffs
import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np

# Handles isomorphism for bonds
bond_match = iso.categorical_edge_match('bond', '')

# Handles isomorphism for atoms with regards to perodic boundary conditions
ads_match = iso.categorical_node_match(['index', 'ads'], [-1, False]) 

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

def bond_symbol(atoms, a1, a2):
    return "{}{}".format(*sorted((atoms[a1].symbol, atoms[a2].symbol)))

def node_symbol(atom, offset):
    return "{}:{}[{},{},{}]".format(atom.symbol, atom.index, offset[0], offset[1], offset[2])

def add_atoms_node(graph, atoms, a1, o1, **kwargs):
    graph.add_node(node_symbol(atoms[a1], o1), index=a1, **kwargs)

def add_atoms_edge(graph, atoms, a1, a2, o1, o2, adsorbate_atoms, **kwargs):
    dist = 2 - (1 if a1 in adsorbate_atoms else 0) - (1 if a2 in adsorbate_atoms else 0)

    graph.add_edge(node_symbol(atoms[a1], o1),
                   node_symbol(atoms[a2], o2),
                   bond=bond_symbol(atoms, a1, a2),
                   index='{}:{}'.format(*sorted([a1, a2])),
                   dist=dist,
                   ads_only=0 if (a1 in adsorbate_atoms and a2 in adsorbate_atoms) else 2,
                   **kwargs)

def compare_chem_envs(chem_envs1, chem_envs2):
    """Compares two sets of chemical environments to see if they are the same
    in chemical identity.  Useful for detecting if two sets of adsorbates are
    in the same configurations.
    
    Args:
        chem_envs1 (list[networkx.Graph]): A list of chemical environments, 
                                           including duplicate adsorbates
        chem_envs2 (list[networkx.Graph]): A list of chemical environments, 
                                           including duplicate adsorbates
    
    Returns:
        bool: Is there a matching graph (site / adsorbate) for each graph?
    """
    # Do they have the same number of adsorbates?
    if len(chem_envs1) != len(chem_envs2):
        return False

    envs_copy = chem_envs2[:] # Make copy of list

    # Check if chem_envs1 matches chem_envs2 by removing from envs_copy
    for env1 in chem_envs1: 
        for env2 in envs_copy:
            if nx.is_isomorphic(env1, env2, edge_match=bond_match):
                # Remove this from envs_copy and move onto next env in chem_envs1
                envs_copy.remove(env2)
                break

    # Everything should have been removed from envs_copy if everything had a match
    if len(envs_copy) > 0:
        return False

    return True

def unique_chem_envs(chem_envs_groups, metadata=None, verbose=True):
    """Given a list of chemical environments, find the unique
    environments and keep track of metadata if required.

    This function exists largely to help with unique site detection
    but its performance will scale badly with extremely large numbers
    of chemical environments to check.  This can be split into parallel
    jobs.

    Args:
        chem_env_groups (list[list[networkx.Graph]]): 
            Chemical environments to compare against each other
        metadata (list[object]): 
            Corresponding metadata to keep with each chemical environment

    Returns:
        list[list[list[networkx.Graph]]]: A list of unique chemical environments 
                                          with their duplicates
        list[list[object]]: A matching list of metadata
    """
    # Error checking, this should never really happen
    if len(chem_envs_groups) == 0:
        return [[],[]]

    # We have metadata to manage
    if metadata is not None:
        if len(chem_envs_groups) != len(metadata):
            raise ValueError("Metadata must be the same length as\
                              the number of chem_envs_groups")
    
    # No metadata to keep track of
    if metadata is None:
        metadata = [None] * len(chem_envs_groups)

    # Keep track of known unique environments
    unique = []

    for index, env in enumerate(chem_envs_groups):
        for index2, (unique_indices, unique_env) in enumerate(unique):
            if verbose:
                print("Checking for uniqueness {:05d}/{:05d} {:05d}/{:05d}".format(index+1, len(chem_envs_groups), index2, len(unique)), end='\r')
            if compare_chem_envs(env, unique_env):
                unique_indices.append(index)
                break
        else: # Was unique
            unique.append(([index], env))
    
    # Zip trick to split into two lists to return
    return zip(*[(env, [metadata[index] for index in indices]) for (indices, env) in unique])

def unique_adsorbates(chem_envs):
    """Removes duplicate adsorbates which occur when perodic
    boundary conditions detect the same adsorbate in two places. Each
    adsorbate graph has its atomic index checked to detect when PBC has 
    created duplicates.

    Args:
        chem_env_groups (list[networkx.Graph]): Adsorbate/environment graphs

    Returns:
        list[networkx.Graph]: The unique adsorbate graphs
    """
    # Keep track of known unique adsorbate graphs
    unique = []
    for env in chem_envs:
        for unique_env in unique:
            if nx.is_isomorphic(env, unique_env, edge_match=bond_match, node_match=ads_match):
                break
        else: # Was unique
            unique.append(env)
    return unique


def process_site(atoms, full, nl, site, radius=3):
    #neighbors = nl.neighbors
    #offsets = nl.displacements
    #neighbors, offsets = nl.get_neighbors()
    full.add_node("X", index=None)
    offset = np.array([0, 0, 0])
    full.add_edge("X",
                  node_symbol(atoms[site[0]], offset),
                  bond="X:{}".format(atoms[site[0]].symbol),
                  ads=0) # Handle first manually
    for last_index, next_index in zip(site[:-1], site[1:]):
        # Error handling needed, .index could be None / -1?
        neighbors, offsets = nl.get_neighbors(last_index)
        #neighbor_index = list(neighbors[last_index]).index(next_index)
        neighbor_index = list(neighbors).index(next_index)
        #offset += offsets[last_index][neighbor_index]
        offset += offsets[neighbor_index]
        #print(offset)
        full.add_edge("X",
                      node_symbol(atoms[next_index], offset),
                      bond="X:{}".format(atoms[next_index].symbol),
                      ads=0)

    site_graph = nx.ego_graph(full, "X", radius=(radius*2)+1, distance="dist")
    site_graph = nx.subgraph(full, list(site_graph.nodes()))
    site_graph = nx.Graph(site_graph)
    full.remove_node("X")
    return site_graph

def process_atoms(atoms, nl, adsorbate_atoms=None, radius=2, grid_n=(2, 2, 0), clean_graph=None):
    """Takes an ase Atoms object and processes it into a full graph as well as
    a list of adsorbate graphs.  This allows for the further post processing
    of the graph to identify chemical environments, chemical identity, and
    site detection.

    Args:
        atoms (ase.Atoms): The atoms object to process
        nl (ase.neighborlist.Neighborlist): A neighborlist built on the atoms object
        adsorbate_atoms (list[int]): The indices of adsorbate atoms
        radius (int): The radius for adsorbate graphs, this is a tunable parameter
        grid_n (tuple[int]): (X,Y,Z) repetitions for PBC.  This should be raised until
                             graphs do not form loops across PBC.  A good starting point
                             for surface DFT is (2,2,0) but this may change if radius is
                             large.

    Returns:
        networkx.Graph: The full graph containing every atom
        list[networkx.Graph]: A list of adsorbate graphs without duplicates from perodic
                              boundary conditions.
    """
    distances = atoms.get_all_distances(mic=True)

    full = nx.Graph()

    # Grid argument determines how many edge repetitions are added.  Scaling for this will be extremely bad
    # (2 * (n - 1) + 1) ^ 2 scaling

    # Add all atoms to graph
    for index, atom in enumerate(atoms):
        for x, y, z in grid_iterator(grid_n):
            add_atoms_node(full, atoms, index, (x, y, z))   

    # Add all edges to graph
    for index, atom in enumerate(atoms):
        for x, y, z in grid_iterator(grid_n):
            neighbors, offsets = nl.get_neighbors(index)
            for neighbor, offset in zip(neighbors, offsets):
                ox, oy, oz = offset
                if not (-grid_n[0] <= ox + x <= grid_n[0]):
                    continue
                if not (-grid_n[1] <= oy + y <= grid_n[1]):
                    continue
                if not (-grid_n[2] <= oz + z <= grid_n[2]):
                    continue
                # This line ensures that only surface adsorbate bonds are accounted for that are less than 2.5 Å
                if distances[index][neighbor] > 2.5 and (bool(index in adsorbate_atoms) ^ bool(neighbor in adsorbate_atoms)):
                    continue
                add_atoms_edge(full, atoms, index, neighbor, (x, y, z), (x + ox, y + oy, z + oz), adsorbate_atoms)

    # Add the case of surface-ads + ads-ads edges to clean graph case here
    if clean_graph:
        edges = [(u, v, d) for u, v,d in full.edges.data() if d["dist"] < 2]    ### Read all the edges, that are between adsorbate and surface (dist<2 condition)
        nodes = [(n, d) for n, d in full.nodes.data() if d["index"] in adsorbate_atoms]    ### take all the nodes that have an adsorbate atoms
        full=nx.Graph(clean_graph)
        full.add_nodes_from(nodes)
        full.add_edges_from(edges)
        
    
    # All adsorbates into single graph, no surface
    ads_nodes = None
    ads_nodes = [node_symbol(atoms[index], (0, 0, 0)) for index in adsorbate_atoms]
    ads_graphs = nx.subgraph(full, ads_nodes)

    # Cut apart graph
    ads_graphs = connected_component_subgraphs(ads_graphs) 

    chem_envs = []
    for ads in ads_graphs:
        initial = list(ads.nodes())[0]
        full_ads = nx.ego_graph(full, initial, radius=0, distance="ads_only")

        new_ads = nx.ego_graph(full, initial, radius=(radius*2)+1, distance="dist")
        new_ads = nx.Graph(nx.subgraph(full, list(new_ads.nodes())))

        for node in full_ads.nodes():
            new_ads.add_node(node, ads=True)
        
        chem_envs.append(new_ads)

    chem_envs = unique_adsorbates(chem_envs)  

    chem_envs.sort(key=lambda x: len(x.edges()))

    return full, chem_envs


if __name__ == "__main__":
    from sys import argv
    from ase.io import read
    import argparse

    parser = argparse.ArgumentParser(
        description='Some command line utilities for atomic graphs')
    parser.add_argument(
        '--radius', '-radius', '-r',
        type=int, default=2,
        help='Sets the graph radius, this can be tuned for different behavior')
    parser.add_argument(
        '--mult', '-mult', '-m',
        type=float, default=1.1,
        help='Sets the radii multiplier for the neighborlist')
    parser.add_argument(
        '--skin', '-skin', '-s',
        type=float, default=0.25,
        help='Sets the skin for the neighborlist, \
              this is a constant offset extending a bond')
    parser.add_argument(
        '--adsorbate-atoms', '-adsorbate-atoms', '-a',
        type=str, default='',
        help='Comma delimited list of elements to be considered adsorbate')
    parser.add_argument(
        '--grid', '-grid', '-g',
        type=str, default='2,2,0',
        help='Grid for PBC as comma delimited list, default is for surface')
    parser.add_argument(
        '--view-adsorbates', '-view-adsorbates', 
        '--view', '-view', '-v',
        action='store_const',
        const=True, default=False,
        help='Displays the adsorbates using matplotlib')
    parser.add_argument(
        '--unique', '-unique', '-u',
        action='store_const',
        const=True, default=False,
        help='Outputs the unique chemical environments')
    parser.add_argument(
        'filenames',
        type=str,
        nargs='+',
        help='Atoms objects to process')
    parser.add_argument(
        '--clean', '-clean','-c', 
        type=str, default=None,
        help='Use clean atoms object')
    args = parser.parse_args()

    radii_multiplier = args.mult
    adsorbate_elements = args.adsorbate_atoms.split(",")

    all_atoms = []
    filenames = []

    for filename in args.filenames:
        try:
            all_atoms.append(read(filename))
            filenames.append(filename)
        except Exception as e:
            print("{} failed to read".format(filename))

    chem_envs = []
    energies = []
    if args.clean:
        clean_atoms = read(args.clean)
        nl = NeighborList(natural_cutoffs(clean_atoms, radii_multiplier), self_interaction=False, 
                          bothways=True, skin=args.skin)
        nl.update(clean_atoms)
        adsorbate_atoms = [atom.index for atom in clean_atoms if atom.symbol in adsorbate_elements]
        if len(adsorbate_atoms):
            raise Exception('Clean atoms should not have adsorbate atoms')
        args.clean, chem_env = process_atoms(clean_atoms, nl, adsorbate_atoms=adsorbate_atoms, 
                                    radius=args.radius, grid_n=[int(grid) for grid in args.grid.split(",")])   ## store the full graph for clean surface, to be used later

    for atoms in all_atoms:
        try:
            energies.append(atoms.get_potential_energy())  ## Attempt to read the OUTCAR here
        except:
            energies.append(float('inf'))

        nl = NeighborList(natural_cutoffs(atoms, radii_multiplier), self_interaction=False, 
                          bothways=True, skin=args.skin)
        nl.update(atoms)
        adsorbate_atoms = [atom.index for atom in atoms if atom.symbol in adsorbate_elements]
        full, chem_env = process_atoms(atoms, nl, adsorbate_atoms=adsorbate_atoms, 
                                       radius=args.radius, grid_n=[int(grid) for grid in args.grid.split(",")],clean_graph=args.clean)  ## get the full graph and the chem_env for all the adsorbates found
        chem_envs.append(chem_env)

        labels = [None]*len(chem_env)
        for index, graph in enumerate(chem_env):
            labels[index] = {node:str(len([edge for edge in full[node] if edge.split(":")[0] not in adsorbate_elements])) for node in graph.nodes()}

    ###### This next condition, finds the unique configs amongst the OUTCARS/ atoms object provided and arranges them according to ascending order of energies
    if args.unique:
        unique, groups = unique_chem_envs(chem_envs, list(zip(energies, filenames)))
        print("Outputting the lowest energy unique configurations")
        groups = [sorted(group) for group in groups]
        for group in sorted(groups):
            if group[0][0] == float('inf'):
                print("{}: Unknown energy, {} duplicates".format(group[0][1], len(group) - 1))
            else:
                print("{}: {} eV, {} duplicates".format(group[0][1], group[0][0], len(group) - 1))
            for duplicate in group[1:]:
                if duplicate[0] == float('inf'):
                    print("-> {}: Unknown energy".format(duplicate[1]))
                else:
                    print("-> {}: {} eV".format(duplicate[1], duplicate[0]))
    if args.view_adsorbates:
        print('Opening graph for viewing')
        #print(chem_env)
        #print(labels)
        draw_atomic_graphs(chem_env, atoms=atoms, labels=labels)