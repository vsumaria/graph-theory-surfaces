==============
Surface Graphs: General Settings and Utilities
==============

The general settings can be divided in the following categories:
 * Graph generation settings
 * Graph analysis settings
 * Graph post-processing for different applications

*************************
Graph generation settings
*************************
These are settings that the user need to keep in mind, while attempting to generate a graph for a given system
 * Neigborlist: To generate the graph, the package requires a neighborlist. We currently use the neighborlist generated using the Atomic Simulation Environment (https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html). For the case of ase-neighborlist, the parameter of scaling the atomic radii and the skin factor are the required inputs. 
  * Atomic Radii Scaling: This number specifies how much the covalent radii of a given atom should be scaled. We have found a modest scaling of '1.1-1.2' to be sufficient for metal slabs. 
  * skin factor: This is an input used as an error factor for all the bonds. The default setting of 0.25 Angstroms is sufficient.
 * Grid parameter: This parameter represents the 
