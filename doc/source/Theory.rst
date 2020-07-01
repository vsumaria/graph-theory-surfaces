==============
Surface Graphs
==============

Surface graphs are made to accomplish two primary tasks.

* Unwrap the full connectivity of the cell into a repeating non-periodic graph to simplify further analysis and make visualization easier.
* Isolate adsorbates and establish the local chemical environment based on a distance criteria.
 * This can be used for further post processing analysis for multidentate adsorbates.
 * This can also be used to compare different calculations for uniqueness of adsorbate configurations.

******************
Full Surface Graph
******************

Initially, a surface graph must be created that encompasses the entire unit cell.  

First, nodes are added for each atom within the cell as well as repetitions in a user specified grid (grid_n) in order to unwrap periodic boundary conditions.  Then edges are added for each bond ("bond") and the bonds are labelled in the following format with the elements in alphabetical order: AB (Example: PtSn). Edges are also assigned two distances for chemical environment analysis.  The distance ("dist") is defined as 0, 1, or 2 based on how many surface atoms are involved as well as an extra distance ("ads_only") which is either 0 or 2.

* If the bond is between ads-ads, both distances are set to 0.
* If the bond is between surface-surface, both distances are set to 2.
* If the bond is between ads-surface, "dist" is set to 1 and "ads_only" is set to 2.

**TODO: Review if "ads_only" is actually needed.**

*********************
Chemical Environments
*********************

The chemical environment is defined as the ego graph with a given radius which is representative of the number of surface coordination shells that should be captured.  A shell radius of 0 will capture only the adsorbate and its direct connections to the surface.  A shell radius of >1 will capture that many coordination shells on the surface.  This is implemented by taking and ego graph from some atom of the adsorbate with a scaled radius.  The graph radius is given as twice the shell radius plus one.  This results in the following behavior.

1. First the entire adsorbate is captured for free since the distance is given as 0.
2. The plus one of the radius will then capture the surface atoms bonded to the adsorbate, given that this is almost never not desired.  This should result in an even radius leftover for capturing shells on the surface.
3. Surface atoms will then be captured by decrementing the radius by 2 each time.  This should always result in an even number of leftover radius.
4. If another adsorbate is found, the leftover radius is decremented by 1, the entire adsorbate is captured for free, then it is decremeneted by 1 again as it captures all the surface atoms.  This also results in an even number of leftover radius.

Then a subgraph of the full graph is taken that contains all of the nodes found in the ego graph.  This effectively ensures that all edges between all found nodes are captured, even if they were not traced out directly in the graph.  This prevents some dangling bonds at the edges of the nodes and gives some additional information that we have found to be useful, otherwise an extra shell is required to find roughly the same results.

Knowing what an appropriate radius is depends on the system and should be tested.  For example, we have found that in FCC metals that a radius of 2 is often sufficient and not much effect is seen by going above 3.
