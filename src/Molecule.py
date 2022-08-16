
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
    
        ##### Miscellaneous Molecular Data #####
        self.amino_acid: bool = len(AMINO_ACID_REGEX.findall(smiles)) != 0 
            
    def createEdges(self) -> "list[Edge]":
        """Create the edges of the software molecule graph using the SMILES code"""

        ##### Algorithm Variables #####
        atom_index: int = 0
        match_index: int = 0
        ring_queue: "dict[str, int]" = {}
        parenth_stack: list[int] = []
        bond: Literal["", "=", "#"] = ""
        edges: "list[Edge]" = []
        
        ##### Algorithm Implementation #####
        for symbol in self.smiles[1:]:

            ##### Atom Symbol Case #####
            if ATOM_REGEX.match(symbol):
                atom_index+=1
                edge_atoms = [self.verticies[match_index], self.verticies[atom_index]]
                new_edge = Edge(edge_atoms, bond)
                edges.append(new_edge)
                match_index = atom_index
                bond = ""

            ##### Bond Symbol Case #####
            if BOND_REGEX.match(symbol):
                bond = symbol               # type: ignore

            ##### Digit Symbol Case #####
            if DIGIT_REGEX.match(symbol):
                if symbol in ring_queue:
                    ring_atom_index = ring_queue.pop(symbol)
                    edge_atoms = [self.verticies[ring_atom_index], self.verticies[atom_index]]
                    new_edge = Edge(edge_atoms, "")
                    edges.append(new_edge)
                else:
                    ring_queue[symbol] = atom_index 

            ##### Parenthesis Symbol Case #####
            if PARENTH_REGEX.match(symbol):
                if symbol == '(':
                    parenth_stack.append(atom_index)
                else:
                    match_index = parenth_stack.pop()
        
        ##### Algorithm Check #####
        assert not parenth_stack
        assert not ring_queue

        ##### Algorithm Results #####
        return edges

        
    def createRings(self):
        """ Determine the number of aromatic and non-aromatic rings using the SMILES code
            Distinguish all atoms as aromatic, non-aromatic, or non-cyclic using the SMILES code
        """

        ########## Algorithm Preparation ##########

        ##### Preparation Variables #####
        ring_index: int = 0
        p_group_counter: int = 0
        parenth_stack: list[int] = [0]
        ring_queue: dict[str, int] = {}
        ring_info: dict[int, list[int]] = {}

        ##### Preparation Implementation #####
        for symbol in self.smiles[1:]:

            ##### Digit Symbol Case #####
            if DIGIT_REGEX.match(symbol):

                if symbol in ring_queue:
                    ring_queue.pop(symbol)

                else:
                    ring_queue[symbol] = ring_index
                    ring_info[ring_index] = [parenth_stack[-1]]
                    ring_index+=1

            ##### Parenthesis Symbol Case #####
            if PARENTH_REGEX.match(symbol):

                if symbol == '(':
                    p_group_counter+=1
                    parenth_stack.append(p_group_counter)
                    for ring_idx in ring_queue.values():
                        ring_info[ring_idx].append(p_group_counter)

                else:
                    closing_p_group = parenth_stack.pop(-1)
                    for ring_idx in ring_queue.values():
                        ring_info[ring_idx] = [p_group for p_group in ring_info[ring_idx] if p_group != closing_p_group]

        ##### Preparation Check #####
        assert not ring_queue
        assert parenth_stack == [0]

        ##### Preparation Results #####
        # print(ring_info)

        ########## Algorithm Implementation ##########

        ##### Algorithm Variables #####
        ring_index = 0
        p_group_counter = 0
        atom_index: int = 0
        ring_stack: list[int] = []
        ring_set: dict[int, list[int]] = {}
        ring_p_groups: set[int] = set()
        ring_atom_indices: list[int] = []

        ##### Algorithm Implementation #####
        for symbol in self.smiles[1:]:

            ##### Atom Symbol Case #####
            if ATOM_REGEX.match(symbol):
                atom_index+=1

                if ring_queue:

                    if parenth_stack[-1] in ring_p_groups:
                        ring_atom_indices.append(atom_index)

                    if parenth_stack[-1] in ring_info[ring_stack[-1]]:
                        ring_set[ring_stack[-1]].append(atom_index)

            ##### Digit Symbol Case #####
            if DIGIT_REGEX.match(symbol):

                if symbol in ring_queue:

                    close_ring_index = ring_queue.pop(symbol)

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
                    ring_atom_indices.append(atom_index)
                    ring_index+=1
                
                ring_stack = list(ring_queue.values())
                ring_p_groups = set(
                    itertools.chain.from_iterable(
                        [p_groups for ring_idx, p_groups in ring_info.items() if ring_idx in ring_queue.values()]
                    )
                )


            ##### Parenthesis Symbol Case #####
            if PARENTH_REGEX.match(symbol):

                if symbol == '(':
                    p_group_counter+=1
                    parenth_stack.append(p_group_counter)

                else:
                    parenth_stack.pop(-1)


        ##### Algorithm Check #####
        assert not ring_queue
        assert parenth_stack == [0]

        ##### Algorithm Results #####
        # print(ring_set)
        # print(ring_atom_indices)

        ########## Algorithm Collection ##########

        ##### Collection Variables #####
        aromatic_ring_count: int = 0
        non_aromatic_ring_count: int = 0

        ##### Collection 1: Ring Counts #####
        for (ring_idx, atom_indices) in ring_set.items():

            if len([v for v in atom_indices if self.verticies[v].ring_type == "aromatic"]) == len(atom_indices):
                aromatic_ring_count+=1
            else:
                non_aromatic_ring_count+=1

        ##### Collection 2: Atom Ring Types #####
        for (ring_idx, atom_indices) in ring_set.items():

            for atom_index in atom_indices:
                if self.verticies[atom_index].ring_type == "non-cyclic":
                    self.verticies[atom_index].ring_type = "non-aromatic"

        ##### Collection Check #####
        assert len(ring_info.keys()) == (aromatic_ring_count + non_aromatic_ring_count)

        ##### Collection Results #####
        return (
            set(ring_atom_indices),  
            {
                "aromatic": aromatic_ring_count, 
                "non-aromatic": non_aromatic_ring_count
            }
        )


    def __str__(self):
        return self.name + " : " + ''.join(self.smiles)

    def __repr__(self):
        return self.name + " : " + ''.join(self.smiles)


