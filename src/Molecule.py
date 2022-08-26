#! 3.10.1
from collections import Counter, defaultdict
from copy import copy
import itertools
from typing import Literal
from Vertex import Vertex
from Edge import Edge
from config import FGSPATH
from constants import ATOM_REGEX, BOND_REGEX, CHARGE_REGEX, DIGIT_REGEX, PARENTH_REGEX, AMINO_ACID_REGEX, SMILES_REGEX, VALENCE_COUNTS, ELECTRON_BOND_COUNTS



class Molecule():
    """ Defines the software molecule space for an organic molecule based on the SMILES code
        through creation of a Simple Undirected Connected Molecular Graph and an identification of ring atoms
    """

    def __init__(self, smiles: str, name: str, type: Literal["mol", "fg"]):
        """ Initialize a molecule obejct from smiles with software molecular graph and identified rings"""

        ##### Input Data #####
        self.smiles: list[str] = [symbol for symbol in SMILES_REGEX.findall(smiles)]
        self.atoms: list[str] = ATOM_REGEX.findall(smiles)
        self.name: str = name
        assert ['[', ']'] not in self.atoms

        ##### Software Molecule Graph (Graph Theory) #####
        self.verticies: "list[Vertex]" = self.createVertices()
        self.order: int = len(self.verticies)
        self.edges: "list[Edge]" = self.createEdges() 
        self.size: int = len(self.edges)
        assert self.order == len(self.atoms)

        ##### Ring Data #####
        self.ring_atoms, self.aromatic_ring_count, self.non_aromatic_ring_count = self.createRings()
        self.total_rings: int = self.aromatic_ring_count + self.non_aromatic_ring_count
        self.total_ring_atoms: int = len(self.ring_atoms)
        self.total_aromatic_atoms: int = len([symbol for symbol in self.atoms if symbol.islower()])
        self.total_non_aromatic_atoms: int = self.total_ring_atoms-self.total_aromatic_atoms
    
        ##### Atom Counts #####
        self.atom_counts: dict[str, int] = Counter([v.symbol for v in self.verticies])
    
        ##### Miscellaneous Molecular Data #####
        self.amino_acid: bool = len(AMINO_ACID_REGEX.findall(smiles)) != 0 

        ##### Functional Groups #####
        self.functional_groups_all, self.functional_groups_exact = self.createFunctionalGroups() if type == "mol" else ({}, {})


    def createVertices(self) -> "list[Vertex]":
        """Create the vertices of the software molecule graph using the SMILES code"""

        ##### Vertex List and Objects #####
        vertices: list[Vertex] = []
        vertex: Vertex = Vertex()
        charge: str = ""

        ##### Atom Symbols Loop #####
        for index, symbol in enumerate(self.atoms):

            ##### Charge Symbol Case #####
            if CHARGE_REGEX.search(symbol):
                charge = CHARGE_REGEX.search(symbol).group()               # type: ignore
                symbol = ''.join([char for char in symbol if charge != char])

            ##### Vertex Object Construction #####
            vertex = Vertex(
                index=index, 
                symbol=symbol.upper() + charge, 
                isAromatic=symbol.islower(), 
                valenceElectronsRequired=VALENCE_COUNTS[symbol.upper()],
                charge=charge
            )

            ##### Reset Charge #####
            charge = ""

            ##### Append Vertex Object #####
            vertices.append(vertex)

        ##### Return Verticies #####
        return vertices
            
    def createEdges(self) -> "list[Edge]":
        """Create the edges and the vertex degrees of the software molecule graph using the SMILES code"""

        ##### Algorithm Variables #####
        atom_index: int = 0
        match_index: int = 0
        edge_index: int = 0 
        ring_queue: "dict[str, int]" = {}
        parenth_stack: list[int] = []
        bond: Literal["", "=", "#"] = ""
        edges: "list[Edge]" = []
        
        ##### Algorithm Implementation #####
        for i,symbol in enumerate(self.smiles[1:]):

            ##### Atom Symbol Case #####
            if ATOM_REGEX.match(symbol):
                atom_index+=1
                edge_atoms = [self.verticies[match_index], self.verticies[atom_index]]
                new_edge = Edge(edge_atoms, bond, edge_index)
                edge_index+=1
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
                    new_edge = Edge(edge_atoms, "", edge_index)
                    edge_index+=1
                    edges.append(new_edge)
                else:
                    ring_queue[symbol] = atom_index 

            ##### Parenthesis Symbol Case #####
            if PARENTH_REGEX.match(symbol):
                if symbol == '(':    
                    # double parenthetical groups will re-append the match index
                    if self.smiles[1:][i-1] == ')':
                        parenth_stack.append(match_index)
                    else:
                        parenth_stack.append(atom_index)
                else:
                    match_index = parenth_stack.pop()

        
        ##### Algorithm Check #####
        assert not parenth_stack
        assert not ring_queue


        ##### Set Vertex Degrees #####
        for vertex in self.verticies:

            ##### R Vertex Degree #####
            if vertex.symbol == 'R':
                vertex.implicit_degree = 0
                vertex.explicit_degree = 1
                vertex.total_degree = 1
                continue

            ##### Core Vertex Degree #####
            total_edges = [edge for edge in edges if vertex.index in edge.indices]
            explicit_valence_electrons = sum([ELECTRON_BOND_COUNTS[edge.bond_type] for edge in total_edges])        
            implicit_valence_electrons = vertex.valence_electrons_required - explicit_valence_electrons            # number of hydrogens
            vertex.explicit_degree = len(total_edges)
            vertex.implicit_degree = implicit_valence_electrons
            vertex.total_degree = vertex.implicit_degree + vertex.explicit_degree

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
        for atom_index in ring_atom_indices:
            if self.verticies[atom_index].ring_type == "non-cyclic":
                self.verticies[atom_index].ring_type = "non-aromatic"

        ##### Collection Check #####
        assert len(ring_info.keys()) == (aromatic_ring_count + non_aromatic_ring_count)

        ##### Collection Results #####
        return (
            set(ring_atom_indices),  
            aromatic_ring_count,
            non_aromatic_ring_count,
        )


    def createFunctionalGroups(self):
        """Wrapper for DFS"""

        all_fgs: list[Molecule] = []
        for (fg_smiles, fg_name) in [                                               
            [y.strip() for y in x.split(' ')]                                       
            for x                                                                   
            in open(FGSPATH.resolve(), "r+").readlines()                            
        ]:

            fg = Molecule(fg_smiles, fg_name, "fg")
            fg_matches = []   

            matched_core_fg_atoms = {
                fg_vertex.index: [
                    mol_vertex for mol_vertex in self.verticies 
                    if mol_vertex.symbol == fg_vertex.symbol and 
                    mol_vertex.total_degree == fg_vertex.total_degree
                ]
                for fg_vertex in fg.verticies if fg_vertex.symbol != 'R'
            }

            for fg_vertex_index, matched_mol_vertices in matched_core_fg_atoms.items():
                fg_vertex = fg.verticies[fg_vertex_index]
                for mol_vertex in matched_mol_vertices:
                    fg_matched_atoms, _, _, _ = self.DFS(fg, fg_vertex, mol_vertex, [], [], True)
                    # if all core atoms in the fg have been given an exact 1:1 pair with a corresponding organic molecule atom
                    # then an fg match by the DFS algorithm has been made
                    if len(fg_matched_atoms) == len([vertex for vertex in fg.verticies if vertex.symbol != 'R']):
                        fg_match = Molecule(fg_smiles, fg_name, "fg")
                        for (fg_atom_index, om_atom_index) in fg_matched_atoms.items():
                            fg_match.verticies[fg_atom_index].index = self.verticies[om_atom_index].index
                            fg_match.verticies[fg_atom_index].ring_type = self.verticies[om_atom_index].ring_type
                        aromatic_tally, non_aromatic_tally = 0,0
                        for vertex in [fg_vertex for fg_vertex in fg_match.verticies if fg_vertex.symbol != 'R']:
                            if vertex.ring_type == "aromatic":
                                aromatic_tally+=1
                            elif vertex.ring_type == "non-aromatic":
                                non_aromatic_tally+=1
                        if aromatic_tally != 0 or non_aromatic_tally != 0:
                            nomenclature = "Aromatic" if aromatic_tally >= non_aromatic_tally else "Cyclic"
                            fg_match.name = nomenclature + fg_match.name
                        fg_matches.append(fg_match)

            repeat_filtered_fg_matches = self.repetitionFilter(fg_matches)

            for fg in repeat_filtered_fg_matches:
                all_fgs.append(fg)  


        all_fgs = self.heirarchyFlter(all_fgs)
        exact_fgs = self.overlapFilter(all_fgs)

        all_fgs_dict = defaultdict(int)
        for fg in all_fgs:
            all_fgs_dict[fg.name] += 1
        exact_fgs_dict = defaultdict(int)
        for fg in exact_fgs:
            exact_fgs_dict[fg.name] += 1

        return (all_fgs_dict, exact_fgs_dict)

            

    def DFS(self, fg, fg_vertex: Vertex, mol_vertex: Vertex, used_mol_edges: "list[int]", used_fg_edges: "list[int]", is_valid: bool):

        matched_indices = {fg_vertex.index: mol_vertex.index}
        fg_core_edges = [edge for edge in fg.edges if fg_vertex.index in edge.indices and not edge.index in used_fg_edges and not 'R' in edge.symbols]
        om_edges = [edge for edge in self.edges if mol_vertex.index in edge.indices and not edge.index in used_mol_edges]

        if fg_vertex.implicit_degree != 0 and mol_vertex.implicit_degree < fg_vertex.implicit_degree:
            return (matched_indices, used_mol_edges, used_fg_edges, False)

        if not fg_core_edges:
            return ({fg_vertex.index: mol_vertex.index}, used_mol_edges, used_fg_edges, True)

        for fg_edge in fg_core_edges:
            fg_complement_vertex = fg_edge.complement_vertex(fg_vertex.index)
            for om_edge in om_edges:
                if om_edge.index not in used_mol_edges:
                    om_corresponding_vertex = om_edge.complement_vertex(mol_vertex.index)
                    if (
                        om_edge == fg_edge 
                        and 
                        om_corresponding_vertex.total_degree == fg_complement_vertex.total_degree
                    ):
                        path = self.DFS(fg, fg_complement_vertex, om_corresponding_vertex, used_mol_edges + [om_edge.index], used_fg_edges + [fg_edge.index], is_valid)
                        if path[3]:
                            matched_path_atoms, matched_mol_path_edges, matched_fg_path_edges, _ = path
                            for matched_fg_atom, matched_mol_atom in matched_path_atoms.items():
                                matched_indices[matched_fg_atom] = matched_mol_atom
                            
                            for om_edge_index in matched_mol_path_edges:
                                used_mol_edges.append(om_edge_index)
                            
                            for fg_edge_index in matched_fg_path_edges:
                                used_fg_edges.append(fg_edge_index)
                            # path_items = [matched_atoms_from_path, mol_edges_used_during_path, fg_edges_used_during_path]
                            # if no path exists, no updates are made and edges are stil available in mol for next fg edge test
                            # for matched_fg_atom, matched_mol_atom in path.items():
                            #     matched_indices[matched_fg_atom] = matched_mol_atom
                            break
            else:
                return (matched_indices, used_mol_edges, used_fg_edges, False)

        return (matched_indices, used_mol_edges, used_fg_edges, True)


    def heirarchyFlter(self, all_fgs):

        eval_indices: set[int] = set()                              # list of fgs already added to be evaluated

        for i, fg in enumerate(all_fgs):
            for j, fg_compare in enumerate(all_fgs):
                if i == j:
                    continue
                if (
                    set([edge for edge in fg.edges if edge.core_type]) == set([edge for edge in fg_compare.edges if edge.core_type]) 
                    and 
                    set([fg_vertex.index for fg_vertex in fg.verticies if 'R' not in fg_vertex.symbol]) == set([fg_vertex.index for fg_vertex in fg_compare.verticies if 'R' not in fg_vertex.symbol])
                ):
                    eval_indices.add(i)
                    eval_indices.add(j)
        
        skip_indices: set[int] = set()

        for i in eval_indices:
            fg: Molecule = all_fgs[i]
            for core_atom in [fg_vertex for fg_vertex in fg.verticies if fg_vertex.symbol != 'R']:
                if (
                    core_atom.explicit_degree == self.verticies[core_atom.index].explicit_degree 
                    and 
                    core_atom.implicit_degree == self.verticies[core_atom.index].implicit_degree
                ):
                    continue
                else:
                    skip_indices.add(i)          # type: ignore
                    break

        return [fg for i, fg in enumerate(all_fgs) if not i in skip_indices]
        

    def overlapFilter(self, all_fgs):
        skip_indices: set[int] = set()
        for i, fg in enumerate(all_fgs):
            for fg_compare in all_fgs:
                if (
                    len([fg_vertex for fg_vertex in fg.verticies if fg_vertex.symbol != 'R']) < 
                    len([fg_vertex for fg_vertex in fg_compare.verticies if fg_vertex.symbol != 'R'])
                ):
                    if set([fg_vertex.index for fg_vertex in fg.verticies if fg_vertex.symbol != 'R']).issubset(set([fg_vertex.index for fg_vertex in fg_compare.verticies if fg_vertex.symbol != 'R'])):
                        skip_indices.add(i)
        return [fg for i, fg in enumerate(all_fgs) if not i in skip_indices]

    def repetitionFilter(self, fg_matches):
        repeat_filtered_fg_matches: list[Molecule] = []
        for fg in fg_matches:
            for already_appeared_fg in repeat_filtered_fg_matches:
                if set([fg_vertex.index for fg_vertex in fg.verticies if not 'R' in fg_vertex.symbol]) == set([fg_vertex.index for fg_vertex in already_appeared_fg.verticies if not 'R' in fg_vertex.symbol]):
                    break
            else:
                repeat_filtered_fg_matches.append(fg)
        return repeat_filtered_fg_matches

    def __str__(self):
        return self.name + " : " + ''.join(self.smiles)

    def __repr__(self):
        return self.name + " : " + ''.join(self.smiles)