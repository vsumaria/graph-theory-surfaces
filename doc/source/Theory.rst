==============
Surface Graphs
==============

Surface graphs are made to accomplish the following tasks:

* Unwrap the full connectivity for a surface type atomic object into a repeating non-periodic graph to simplify further analysis and make visualization easier.
* These graph representations can then identify unique locations on the graph as a function of the chemical environment of the given site. This can then be used:
 * To identify unique starting positions for an adsorbate on a given graph.
 * To identify unique adsorbate configurations amongst a set of configurations.
These two capabilities can be synergistically utilized to enumerate high coverage configurations, as well as, for finding unique configurations post - simulations.  


******************
Full Surface Graph
******************

Initially, a surface graph must be created that encompasses the entire unit cell.  

First, nodes are added for each atom within the cell as well as repetitions in a user specified grid (grid_n) in order to unwrap periodic boundary conditions. 'grid_n' dictates how many periodic repetitions are desired, though we have found that a large grid is not favorable from a computing stand point, and for cell sizes greater or equal to a 3x3, 1-2 repetitions are sufficient. Then edges are added for each bond ("bond") and the bonds are labelled in the following format with the elements in alphabetical order: AB (Example: PtSn). Edges are also assigned two distances for chemical environment analysis.  The distance ("dist") is defined as 0, 1, or 2 based on how many surface atoms are involved as well as an extra distance ("ads_only") which is either 0 or 2.

* If the bond is between ads-ads, both distances are set to 0.
* If the bond is between surface-surface, both distances are set to 2.
* If the bond is between ads-surface, "dist" is set to 1 and "ads_only" is set to 2.

Graphs offer the flexibility to add other important properties to the nodes and the edges as desired by the user. For example, the distance between the surface atoms can be a useful edge property to store, or the coordination number of the surface site can be another important node property to store. Some interesting node and edge properties are listed as follows:
 
* Node Properties:
   * Coordination number
   * D-band center
   * Electronegativity
   * number of valence electrons

* Edge Properties:
   * Atom-atom bond distance
   * Type of bond: surface-surface, surface-ads, ads-ads, H--bond (some of this is captured in the "dist" parameter above)
   * bond energy of a given bond: this can be used to store entities like ads_energy between ads and surface atoms, vacancy formation energy between surface-surface atoms.

**TODO: Review if "ads_only" is actually needed.**

**TODO: Tristan explain how the radius of the graph makes use of this parameter, and how 1.5 captures, the 1st shell surface atoms and then only the ads attached to them**

*********************
Chemical Environments
*********************

The chemical environment is defined as the ego graph with a given radius which is representative of the number of surface coordination shells that should be captured.  A shell radius of 0 will capture only the adsorbate and its direct connections to the surface.  A shell radius of >1 will capture that many coordination shells on the surface.  This is implemented by taking and ego graph from some atom of the adsorbate with a scaled radius.  The graph radius is given as twice the shell radius plus one (r = 2*r_s +1).  This results in the following behavior.

1. First the entire adsorbate is captured for free since the distance is given as 0.
2. The plus one of the radius will then capture the surface atoms bonded to the adsorbate.  This should result in an even radius leftover for capturing shells on the surface (*Can you explain this better?*). Therfore a radius of '0' will provide the adsorbate atoms + the surface atoms surrounding the adsorbate.
3. Surface atoms will then be captured by decrementing the radius by 2 each time.  This should always result in an even number of leftover radius. ** Tristan this left ovr radius is confusing, can you try to explain this a little better, I think the idea is that, for a given user input radius 'r_s', the graph will attempt to find everything in 2*r_s+1, according to the dist parameter right?*
4. If another adsorbate is found, the leftover radius is decremented by 1, the entire adsorbate is captured for free, then it is decremeneted by 1 again as it captures all the surface atoms.  This also results in an even number of leftover radius.  ** This might be the perfect place to explain radius = 1.5 **

Then a subgraph of the full graph is taken that contains all of the nodes found in the ego graph.  This effectively ensures that all edges between all found nodes are captured, even if they were not traced out directly in the graph.  This prevents some dangling bonds at the edges of the nodes and gives some additional information that we have found to be useful, otherwise an extra shell is required to find roughly the same results.

Knowing what an appropriate radius is depends on the system and should be tested.  For example, we have found that in FCC metals that a radius of 2 is often sufficient and not much effect is seen by going above 3. Though this choice can be very system dependent.
