Theoretical Overview
====================

This project makes usage of the molecular graph encoding scheme of the SMILES code to execute a sub-graph 
search algorithm for functional group instances inside of an organic molecule. 

Graphical Theory
----------------

In graph theory, a set of vertices and edges mathematically construct the circles and lines of a graph visualized in space [1]_.
There are various types of graphs, and they are dependent upon the way in which edges connect vertices together, as well the type of data stored in the graph vertices and edges. 
For the SMILES code, the type of encoded graph is a Simple Undirected Connected Molecular graph. The molecular aspect indicates that the data of the vertices 
and edges of this graph are molecular in nature. The simple, undirected, and connected aspect indicates that the edges and vertices follow a specific uniform pattern 
of connection throughout the SMILES derived graph (View the reference for details about this classification [1]_). 
SMILES derived graphs may be functional groups or organic molecules as shown below:


.. figure:: _static/Molecular_Graphs.png

    Organic Molecule and Functional Group Molecular Graphs Example



Functional Groups Theory
------------------------

Functional groups are a collection of sub-structure patterns found in organic molecules that have an influence
on the chemical behaviors of a molecule when exhibited. For example, combinations of functional groups in an organic 
molecule may have energy-efficiency related effects. This type of correllation can be beneficial for opto-electronic 
materials developers in pursuit of reducing energy consumption of opto-electronic devices. 

Functional group sub-structures are charecterized by their `core` structure and their `R` structure. 
The set of edges and vertices which involve only `core` atoms may be called its `core` structure.
The set of edges and vertices which involve `R` atoms may be called its `R` structure. The `R` structure is dependent
upon the `core` structure for every functional group.


Functional groups may be hierarchically related to one another when two or more functional groups have equivalent `core` structures, each with different `R` structures. 
A common example of this behavior is shown for the Amine type functional groups:

.. figure:: _static/Hierarchical_Functional_Groups.png

    Functional Group Hierarchy of Primary, Secondary, and Tertiary Amines

Functional groups may overlap with one another when the `core` structure of one functional group is entirely encapsulated inside the `core` structure of another functional group.
A common exmaple of this behavior is shown for the Ester and Ether functional groups:


.. figure:: _static/Overlap_Functional_Groups.png

    Functional Group Overlap of an Ester with an Ether 

`R` atoms may be generalized as any bond which fills a single valence electron of a `core` atom.
For hierarchically related functional groups, the difference between an `R` atom as a hydrogen or an `R` atom as a non-hydrogen becomes significant.

Functional group sub-structures may be catergorized by their exsistence inside of ring structures. They may be given
the more specific nomenclature of `Aromatic`, `Non Aromatic`, and `Non Cyclic` based on the number of `core` atoms 
which appear inside of aromatic or non aromatic ring structures of the organic molecule. 


.. rubric:: Footnotes
.. [1] Introduction to Graph Theory Fourth Edition Robin J. Wilson 