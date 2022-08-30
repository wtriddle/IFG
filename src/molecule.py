"""A class for decoding chemical information from the Simplified Molecular Input Line Entry System (SMILES)"""

from collections import Counter, defaultdict
import itertools
from typing import Literal
from vertex import Vertex
from edge import Edge
from config import FUNCTIONAL_GROUPS_PATH
from constants import ATOM_REGEX, BOND_REGEX, CHARGE_REGEX, DIGIT_REGEX, PARENTH_REGEX, AMINO_ACID_REGEX, SMILES_REGEX, VALENCE_COUNTS, ELECTRON_BOND_COUNTS

class Molecule():
    """ Defines a Simple Connected Undirected Graph for a molecule using a SMILES code
        Determines the number of aromatic and non-aromatic rings of a molecule using a SMILES code
        Determines the number of unique ring-classified instances of functional groups using the software molecule graph format
    """

    def __init__(self, smiles: str, name: str, type: Literal["mol", "fg"]):
        """ Generates a software molecule graph for any smiles code. 
            Generates ring and functional group data for organic molecules
        """

        ##### Input Data #####
        self.smiles: list[str] = [symbol for symbol in SMILES_REGEX.findall(smiles)]
        self.atoms: list[str] = ATOM_REGEX.findall(smiles)
        self.name: str = name
        assert ['[', ']'] not in self.atoms

        ##### Software Molecule Graph (Graph Theory) #####
        self.vertices: "list[Vertex]" = self.createVertices()
        self.order: int = len(self.vertices)
        self.edges: "list[Edge]" = self.createEdges()
        self.size: int = len(self.edges)
        assert self.order == len(self.atoms)

        ##### Ring Data #####
        self.ring_atoms: set[int]
        self.aromatic_ring_count: int
        self.non_aromatic_ring_count: int
        self.ring_atoms, self.aromatic_ring_count, self.non_aromatic_ring_count = self.createRings()
        self.total_ring_count: int = self.aromatic_ring_count + self.non_aromatic_ring_count
        self.total_ring_atom_count: int = len(self.ring_atoms)
        self.total_aromatic_atoms: int = len([symbol for symbol in self.atoms if symbol.islower()])
        self.total_non_aromatic_atoms: int = self.total_ring_atom_count-self.total_aromatic_atoms
    
        ##### Atom Counts #####
        self.atom_counts: dict[str, int] = Counter([v.symbol for v in self.vertices])
    
        ##### Miscellaneous Molecular Data #####
        self.amino_acid: bool = len(AMINO_ACID_REGEX.findall(smiles)) != 0 

        ##### Functional Groups #####
        self.functional_groups_all: dict[str, int]
        self.functional_groups_exact: dict[str, int]
        self.functional_groups_all, self.functional_groups_exact = self.createFunctionalGroups() if type == "mol" else ({}, {})


    def createVertices(self) -> "list[Vertex]":
        """Create the vertices of a software molecule graph using the SMILES code"""

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
                is_aromatic=symbol.islower(), 
                valence_electrons_required=VALENCE_COUNTS[symbol.upper()],
                charge=charge
            )

            ##### Reset Charge #####
            charge = ""

            ##### Append Vertex Object #####
            vertices.append(vertex)

        ##### Return vertices #####
        return vertices
            
    def createEdges(self) -> "list[Edge]":
        """Create the edges and the vertex degrees of a software molecule graph using the SMILES code"""

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
                edge_atoms = [self.vertices[match_index], self.vertices[atom_index]]
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
                    edge_atoms = [self.vertices[ring_atom_index], self.vertices[atom_index]]
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
        assert not bond


        ##### Set Vertex Degrees #####
        for vertex in self.vertices:

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

            if len([v for v in atom_indices if self.vertices[v].ring_type == "aromatic"]) == len(atom_indices):
                aromatic_ring_count+=1
            else:
                non_aromatic_ring_count+=1

        ##### Collection 2: Atom Ring Types #####
        for atom_index in ring_atom_indices:
            if self.vertices[atom_index].ring_type == "non-cyclic":
                self.vertices[atom_index].ring_type = "non-aromatic"

        ##### Collection Check #####
        assert len(ring_info.keys()) == (aromatic_ring_count + non_aromatic_ring_count)

        ##### Collection Results #####
        return (
            set(ring_atom_indices),  
            aromatic_ring_count,
            non_aromatic_ring_count,
        )


    def createFunctionalGroups(self):
        """Determine the number of unique ring classified functional groups from a list of identifiable functional groups using the software molecule graph format"""

        ##### All Functional Group Matches #####
        all_fgs: list[Molecule] = []

        ##### Functional Group Loop #####
        for (fg_smiles, fg_name) in [                                               
            [y.strip() for y in x.split(' ')]                                       
            for x                                                                   
            in open(FUNCTIONAL_GROUPS_PATH.resolve(), "r+").readlines()                            
        ]:

            ##### Functional Group Graph Template #####
            fg: Molecule = Molecule(fg_smiles, fg_name, "fg")

            ##### Functional Group Matches #####
            fg_matches: list[Molecule] = []   

            ##### Functional Group Mol Vertex Start Locations #####
            matched_core_fg_atoms: dict[int, list[Vertex]] = {
                fg_vertex.index: [
                    mol_vertex for mol_vertex in self.vertices 
                    if mol_vertex.symbol == fg_vertex.symbol and 
                    mol_vertex.total_degree == fg_vertex.total_degree
                ]
                for fg_vertex in fg.vertices if fg_vertex.symbol != 'R'
            }

            ##### Functional Group Mol Vertex Start Locations Loop #####
            for fg_vertex_index, matched_mol_vertices in matched_core_fg_atoms.items():

                ##### Functional Group Start Vertex #####
                fg_vertex: Vertex = fg.vertices[fg_vertex_index]

                ##### Molecule Start Vertex Locations Loop #####
                for mol_vertex in matched_mol_vertices:

                    ##### Functional Group DFS Match Algorithm #####
                    fg_matched_atoms, _, _ = self.DFS(fg, fg_vertex, mol_vertex, [], [])

                    ##### Functional Group Match Case #####
                    if len(fg_matched_atoms) == len([vertex for vertex in fg.vertices if vertex.symbol != 'R']):

                        ##### Functional Group Extraction #####
                        fg_match: Molecule = Molecule(fg_smiles, fg_name, "fg")
                        for (fg_atom_index, om_atom_index) in fg_matched_atoms.items():
                            fg_match.vertices[fg_atom_index].index = self.vertices[om_atom_index].index
                            fg_match.vertices[fg_atom_index].ring_type = self.vertices[om_atom_index].ring_type

                        ##### Ring Classification #####
                        aromatic_tally: int = len([fg_vertex for fg_vertex in fg_match.vertices if fg_vertex.symbol != 'R' and fg_vertex.ring_type == "aromatic"])
                        non_aromatic_tally: int = len([fg_vertex for fg_vertex in fg_match.vertices if fg_vertex.symbol != 'R' and fg_vertex.ring_type == "non-aromatic"])
                        if aromatic_tally != 0 or non_aromatic_tally != 0:
                            nomenclature: str = "Aromatic" if aromatic_tally >= non_aromatic_tally else "Cyclic"
                            fg_match.name = nomenclature + fg_match.name

                        ##### Match Add #####                        
                        fg_matches.append(fg_match)

            ##### Repeat Functional Group Match Filter #####
            repeat_filtered_fg_matches: list[Molecule] = self.repetitionFilter(fg_matches)

            ##### All Functional Group Matches Add #####
            for fg in repeat_filtered_fg_matches:
                all_fgs.append(fg)  

        ##### Hierarchical Functional Group Filter #####
        all_fgs: list[Molecule] = self.hierarchyFilter(all_fgs)

        ##### Overlapping Functional Group Filter #####
        exact_fgs: list[Molecule] = self.overlapFilter(all_fgs)

        ##### All Functional Group Counts #####
        all_fgs_dict = defaultdict(int)
        for fg in all_fgs:
            all_fgs_dict[fg.name] += 1

        ##### Exact Functional Group Counts #####
        exact_fgs_dict = defaultdict(int)
        for fg in exact_fgs:
            exact_fgs_dict[fg.name] += 1

        ##### Algorithm Results #####
        return (all_fgs_dict, exact_fgs_dict)

    def DFS(self, fg: "Molecule", fg_vertex: Vertex, mol_vertex: Vertex, used_mol_edges: "list[int]", used_fg_edges: "list[int]"):
        """Search an organic molecule software graph for the presence of a functional group sub-graph structure"""

        ##### New Atom-Pair Backtrack Variable #####
        matched_indices = {fg_vertex.index: mol_vertex.index}

        ##### Edge Sets #####
        fg_core_edges = [edge for edge in fg.edges if fg_vertex.index in edge.indices and not edge.index in used_fg_edges and not 'R' in edge.symbols]
        om_edges = [edge for edge in self.edges if mol_vertex.index in edge.indices and not edge.index in used_mol_edges]

        ##### Implicit Degree Validation #####
        if fg_vertex.implicit_degree != 0 and mol_vertex.implicit_degree < fg_vertex.implicit_degree:
            return ({}, [], [])

        ##### Functional Group End Graph Boundary Case #####
        if not fg_core_edges:
            return ({fg_vertex.index: mol_vertex.index}, used_mol_edges, used_fg_edges)

        ##### Functional Group Core Edge Set Searching #####
        for fg_edge in fg_core_edges:

            ##### Complementary Functional Group Vertex #####
            fg_complement_vertex = fg_edge.complement_vertex(fg_vertex.index)

            ##### Organic Molecule Edge Set Match Attempts #####
            for om_edge in om_edges:

                ##### Unused Organic Molecule Edges #####
                if om_edge.index not in used_mol_edges:

                    ##### Complementary Molecule Vertex #####
                    om_corresponding_vertex = om_edge.complement_vertex(mol_vertex.index)

                    ##### Edge Structure & Complementary Vertex Degree Validation #####
                    if (
                        om_edge == fg_edge 
                        and 
                        om_corresponding_vertex.total_degree == fg_complement_vertex.total_degree
                    ):

                        ##### DFS Recursion #####
                        path = self.DFS(fg, fg_complement_vertex, om_corresponding_vertex, used_mol_edges + [om_edge.index], used_fg_edges + [fg_edge.index])

                        ##### Backtrack Collection #####
                        if all(path):

                            ##### Backtrack Unpacking #####
                            matched_path_atoms, matched_mol_path_edges, matched_fg_path_edges = path

                            ##### Atom Unpacking #####
                            for matched_fg_atom, matched_mol_atom in matched_path_atoms.items():
                                matched_indices[matched_fg_atom] = matched_mol_atom
                            
                            ##### Molecule Edge Unpacking #####
                            for om_edge_index in matched_mol_path_edges:
                                used_mol_edges.append(om_edge_index)
                            
                            ##### Functional Group Edge Unpacking #####
                            for fg_edge_index in matched_fg_path_edges:
                                used_fg_edges.append(fg_edge_index)
                            
                            ##### Satisfied Functional Group Edge #####
                            break

            ##### Unsatisfied Functional Group Edge #####
            else:
                return ({}, [], [])

        ##### All Functional Group Core Edges Satisfied #####
        return (matched_indices, used_mol_edges, used_fg_edges)

    def hierarchyFilter(self, all_fgs) -> "list[Molecule]":
        """Identify and filter hierarchically related functional group matches"""

        ##### Matches List Evaluation Indices #####
        eval_indices: set[int] = set()

        ##### Hierarchical Functional Group Identification #####
        for i, fg in enumerate(all_fgs):
            for j, fg_compare in enumerate(all_fgs):
                if i == j:
                    continue
                if (
                    set([edge for edge in fg.edges if edge.core_type]) == set([edge for edge in fg_compare.edges if edge.core_type]) 
                    and 
                    set([fg_vertex.index for fg_vertex in fg.vertices if 'R' not in fg_vertex.symbol]) == set([fg_vertex.index for fg_vertex in fg_compare.vertices if 'R' not in fg_vertex.symbol])
                ):
                    eval_indices.add(i)
                    eval_indices.add(j)
        
        ##### Indices To-Be Skipped From Matches List #####
        skip_indices: set[int] = set()

        ##### Hierarchical Accuracy Selection #####
        for i in eval_indices:
            fg: Molecule = all_fgs[i]
            for core_atom in [fg_vertex for fg_vertex in fg.vertices if fg_vertex.symbol != 'R']:
                if (
                    core_atom.explicit_degree == self.vertices[core_atom.index].explicit_degree 
                    and 
                    core_atom.implicit_degree == self.vertices[core_atom.index].implicit_degree
                ):
                    continue
                else:
                    skip_indices.add(i)
                    break

        ##### Apply Skips For Accurate Results #####
        return [fg for i, fg in enumerate(all_fgs) if not i in skip_indices]
        
    def overlapFilter(self, all_fgs) -> "list[Molecule]":
        """Identify and filter overlapping functional group matches"""

        ##### Indices To-Be Skipped From Matches List #####
        skip_indices: set[int] = set()

        ##### Overlapping Functional Group Identification and Accuracy Selection #####
        for i, fg in enumerate(all_fgs):
            for fg_compare in all_fgs:
                if (
                    len([fg_vertex for fg_vertex in fg.vertices if fg_vertex.symbol != 'R']) < 
                    len([fg_vertex for fg_vertex in fg_compare.vertices if fg_vertex.symbol != 'R'])
                ):
                    if set([fg_vertex.index for fg_vertex in fg.vertices if fg_vertex.symbol != 'R']).issubset(set([fg_vertex.index for fg_vertex in fg_compare.vertices if fg_vertex.symbol != 'R'])):
                        skip_indices.add(i)

        ##### Apply Skips For Accurate Results #####
        return [fg for i, fg in enumerate(all_fgs) if not i in skip_indices]

    def repetitionFilter(self, fg_matches) -> "list[Molecule]":
        """Identify and filter repeated functional group matches"""

        ##### Repeat Filteres List Of Matches #####
        repeat_filtered_fg_matches: list[Molecule] = []

        ##### Repeat Identification and Removal #####
        for fg in fg_matches:
            for already_appeared_fg in repeat_filtered_fg_matches:
                if set([fg_vertex.index for fg_vertex in fg.vertices if not 'R' in fg_vertex.symbol]) == set([fg_vertex.index for fg_vertex in already_appeared_fg.vertices if not 'R' in fg_vertex.symbol]):
                    break
            else:
                repeat_filtered_fg_matches.append(fg)

        ##### Return Filtered Set For Accurate Results #####
        return repeat_filtered_fg_matches

    def __str__(self):
        """String Representation of a Molecule"""
        return ''.join(self.smiles)

    def __repr__(self):
        """General Representation of a Molecule (Same as String Representation)"""
        return ''.join(self.smiles)