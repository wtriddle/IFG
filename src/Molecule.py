"""
Molecule attribute notes:
    smiles is a symbol list of every smiles symbol minus hydrogen symbols
    atoms is a symbol list of every atom symbol including double letter atoms in the smiles
    name is the name of the molecule of interest
    verticies are the vertex objects in the graph containing atomic data
    order is the number of verticies in the graph (i.e. the number of atoms in the molecule)
    edges are the vertex to vertex connection in the graph (i.e. the undirected lines connecting verticies)
    size are the number of edges in the graph
    ring_atoms are the set of atom indices which belong to rings
    ring_counts is the total number of aromatic and non-aromtic rings in the molecule
    total_rings is the combined number of total rings (aromatic and non-aromtic)
    total_ring_atoms is the number of atoms which are apart of a ring
    total_aromatic_atoms are the total number of aromatic atoms
    total_non_aromatic_atoms are the total number of non_aromatic_atoms
    atom_counts are the total number of each type of atom, double letter atom inclusive
    amino_acid is whether or not an amino acid is present in the molecule


Algorithm Notes:

    Ring Contruction Notes:
    Due to the complexity of non-aromatic ring polycyclism described in the SMILES code,
    the algorithm cannot identify the exact atoms of each ring index individually for non-aromatic atoms
    Instead, it is only capable of identifying the specific atoms of aromtic rings individually,
    and can only identify that atoms are within a non-aromtic ring rather than the exact ring index
    they are apart of.
    This is due to the highly complex expression of polycyclism inside the SMILES code,
    and reaches outside of the bounds of the necessary implementation of the algorithm.
    (i.e all information is still derivable despite this level of accuracy not being achieved in the algorithm)



Molecule Class Notes:
    -Single core atom FGs will be identifiable in IFG
    -Vertex ring_type attribute will identify ring_types of FGs
    -ring_counts and ring_type attributes identified using createRings() algorithm
    -Edges description of software molecule graph refine IFG graph search algorithm to properly identify all FGs (replaces bondData for more refined understanding and conceptual foundation)
    -SMILES code list of symbols minus explict hydrogens and support Double lettered atoms with usage of ATOM_REGEX and SMILES_REGEX 

"""


from collections import Counter
import itertools
from typing import Literal
from Vertex import Vertex
from Edge import Edge
from constants import ATOM_REGEX, BOND_REGEX, DIGIT_REGEX, PARENTH_REGEX, AMINO_ACID_REGEX, SMILES_REGEX



class Molecule():
    """ Defines the software molecule space for an organic molecule based on the SMILES code
        through creation of a Simple Undirected Connected Molecular Graph and an identification of ring atoms
    """

    def __init__(self, smiles: str, name: str):
        """ Initialize a molecule obejct from smiles with software molecular graph and identified rings"""

        ##### Input Data #####
        self.smiles: list[str] = [symbol for symbol in SMILES_REGEX.findall(smiles) if 'H' not in symbol] # remove all hydrogens
        self.atoms: list[str] = ATOM_REGEX.findall(''.join(self.smiles))
        self.name: str = name
        assert ['+', '-', '[', ']', 'H'] not in self.atoms

        ##### Software Molecule Graph (Graph Theory) #####
        self.verticies: "list[Vertex]" = [Vertex(index, symbol.upper(), symbol.islower()) for index, symbol in enumerate(self.atoms)]
        self.order: int = len(self.verticies)
        self.edges: "list[Edge]" = self.createEdges() 
        self.size: int = len(self.edges)
        assert self.order == len(self.atoms)

        ##### Ring Data #####
        self.ring_atoms, self.ring_counts = self.createRings()
        self.total_rings = self.ring_counts["aromatic"] + self.ring_counts["non-aromatic"]
        self.total_ring_atoms = len(self.ring_atoms)
        self.total_aromatic_atoms = len([symbol for symbol in self.atoms if symbol.islower()])
        self.total_non_aromatic_atoms = self.total_ring_atoms-self.total_aromatic_atoms
    
        ##### Atom Counts #####
        self.atom_counts: dict[str, int] = Counter([v.symbol for v in self.verticies])
    
        #### Miscellaneous Molecular Data ####
        self.amino_acid: bool = len(AMINO_ACID_REGEX.findall(smiles)) != 0 
            
    def createEdges(self):
        """Create the edges of the software molecule graph using the SMILES code"""

        atom_index: int = 0
        match_index: int = 0
        ring_queue: "dict[str, int]" = {}      # actual number involved in the ring as well as the atom_index => {ring_number: atom_index}
        parenth_stack: list[int] = []
        bond: Literal["", "=", "#"] = ""
        edges: "list[Edge]" = []
        

        for symbol in self.smiles[1:]:

            # match atom case
            if ATOM_REGEX.match(symbol):

                atom_index+=1
                edge_atoms = [
                    self.verticies[match_index], 
                    self.verticies[atom_index]
                ]
                new_edge = Edge(edge_atoms, bond)
                edges.append(new_edge)
                match_index = atom_index
                bond = ""

                continue

            # match number case
            if DIGIT_REGEX.match(symbol):

                if symbol in ring_queue:
                    corresponding_index = ring_queue.pop(symbol)
                    edge_atoms = [
                        self.verticies[corresponding_index], 
                        self.verticies[atom_index]
                    ]
                    new_edge = Edge(edge_atoms, "")
                    edges.append(new_edge)
                else:
                    ring_queue[symbol] = atom_index 
                
                continue
            

            # match bond case
            if BOND_REGEX.match(symbol):

                bond = symbol               # type: ignore

                continue


            # match parenth case
            if PARENTH_REGEX.match(symbol):

                if symbol == '(':
                    parenth_stack.append(atom_index)
                else:
                    match_index = parenth_stack.pop()
                
                continue
        
        # assert that the parenthesis stack and the ring_queue are empty
        assert not parenth_stack
        assert not ring_queue

        return edges

        
    def createRings(self):
        """ Determine the ring counts and ring types of each atom in the molecule
            ring_stack[-1] === index of most recently opened ring
            p_stack[-1] === index of current p_group during iteration

            Algorithm Steps:
                Step 1: Define per ring index which p_groups are allowable for each ring
                Step 2: Identify that atoms of each ring by index and identify all atoms apart of rings
                Step 3: Identify the number of aromatic and non-aromatic rings, and label atoms as aromatic, non-aromatic and non-cyclic
        """

        # Step 1: Define per ring index which p_groups are allowable for each ring
        ring_index: int = 0                                     # counter for the ring index in the SMILES
        p_stack: list[int] = [0]                                # parenthetical group stack order (always has 0 as root level)
        p_group_counter: int = 0                                # counter for the p_group and p_stack data structure
        ring_queue: dict[str, int] = {}                         # ring_number: ring index pairing for open rings
        ring_info: dict[int, list[int]] = {}                    # ring_index: [p_groups]

        # fill ring_info (set up our NPR and PR defined rings)
        # Generates per ring index the allowable p_groups where atoms apart of their ring may appear
        for symbol in self.smiles[1:]:

            if PARENTH_REGEX.match(symbol):

                if symbol == '(':
                    p_group_counter+=1
                    p_stack.append(p_group_counter)
                    for ring_idx in ring_queue.values():
                        ring_info[ring_idx].append(p_group_counter)

                else:
                    closing_p_group = p_stack.pop(-1)
                    for ring_idx in ring_queue.values():
                        ring_info[ring_idx] = [p_group for p_group in ring_info[ring_idx] if p_group != closing_p_group]

            if DIGIT_REGEX.match(symbol):

                if symbol in ring_queue:
                    ring_queue.pop(symbol)

                else:
                    ring_queue[symbol] = ring_index
                    ring_info[ring_index] = [p_stack[-1]]
                    ring_index+=1

        # ring queue should be empty and p_stack should end on p_group 0
        assert not ring_queue
        assert p_stack == [0]
        p_group_counter = 0
        ring_index = 0


        # Step 2: Determine the aromatic ring indicies and find all atoms apart of a ring
        atom_index: int = 0                                     # current atomic index in the SMILES
        ring_stack: list[int] = []                              # order of open rings in the iteration 2
        ring_set: dict[int, list[int]] = {}                     # ring_index to atom index list of atoms apart of ring index pairing (used for aromatic ring identification)
        ring_p_groups: set[int] = set()                         # set of p_groups for open rings
        total_ring_indices: list[int] = []                      # list of all atoms apart of rings

        for symbol in self.smiles[1:]:

            if PARENTH_REGEX.match(symbol):

                if symbol == '(':
                    p_group_counter+=1
                    p_stack.append(p_group_counter)

                else:
                    p_stack.pop(-1)

            if DIGIT_REGEX.match(symbol):

                if symbol in ring_queue:

                    close_ring_index = ring_queue.pop(symbol)

                    # Polycylic atom when closing a ring that is consecutive with another ring on the same p_group 
                    if ring_queue:
                        prev_ring_index = ring_stack[ring_stack.index(close_ring_index)-1]
                        p_end_group = ring_info[close_ring_index][-1]
                        if p_end_group in ring_info[prev_ring_index]:
                            ring_set[prev_ring_index].append(atom_index)

                    if not atom_index in ring_set[close_ring_index]:
                        ring_set[close_ring_index].append(atom_index)
                else:
                    ring_queue[symbol] = ring_index
                    ring_set[ring_index] = [atom_index]
                    total_ring_indices.append(atom_index)
                    ring_index+=1
                
                # Update the ring stack and ring_p_groups when a new ring is viewed
                ring_stack = list(ring_queue.values())
                ring_p_groups = set(
                    itertools.chain.from_iterable([p_groups for ring_idx, p_groups in ring_info.items() if ring_idx in ring_queue.values()])
                )

            if ATOM_REGEX.match(symbol):
                atom_index+=1

                # if rings are open
                if ring_queue:

                    # track all atoms apart of a ring
                    if p_stack[-1] in ring_p_groups:
                        total_ring_indices.append(atom_index)

                    # add atom to MROR if atom is in valid p_group of MROR
                    if p_stack[-1] in ring_info[ring_stack[-1]]:
                        ring_set[ring_stack[-1]].append(atom_index)

        # ring_queue should be empty and p_stack should be on group 0
        assert not ring_queue
        assert p_stack == [0]


        # Step 3: Use data from step 2 to identify non-aromatic atoms and identify the number of aromatic rings
        aromatic_ring_count: int = 0
        non_aromatic_ring_count: int = 0

        # Data 1: Determine the aromatic and non-aromatic rings present in the molecule
        for (ring_idx, atom_indices) in ring_set.items():

            # all atoms are aromatic in ring means aromatic ring, otherwise nonaromatic ring
            if len([v for v in atom_indices if self.verticies[v].ring_type == "aromatic"]) == len(atom_indices):
                aromatic_ring_count+=1
            else:
                non_aromatic_ring_count+=1

        # Data 2: Determine the aromatic and non-aromatic atoms present in the molecule
        for (ring_idx, atom_indices) in ring_set.items():

            # all aromatic atoms already labeled, if an atom appears in a ring list, then it must be non-aromatic
            for atom_index in atom_indices:
                if self.verticies[atom_index].ring_type == "non-cyclic":
                    self.verticies[atom_index].ring_type = "non-aromatic"

        return (
            set(total_ring_indices),  
            {"aromatic": aromatic_ring_count, "non-aromatic": non_aromatic_ring_count}
        )


    def __str__(self):
        return self.name + " : " + ''.join(self.smiles)

    def __repr__(self):
        return self.name + " : " + ''.join(self.smiles)

