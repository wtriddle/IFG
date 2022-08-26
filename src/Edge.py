"""An Edge is representative of a bond in the molecule, and appears as the edges of a simple undirected connected graph (graph theory)

Edge bond_type is the type of bond between two verticies (i.e atoms) in the molecular graph
Edge verticies are the two atom symbols which are connected by the edge
Edge indices are the two atom index values which identify the two atoms that are connected by the edge (local to the graph)
Edge core_type is the assertion of both verticies being non-R atoms (the negation of core_type is called an R_type which has an R atom involved in the edge connection)
Edge index is the index of the edge in the graph
"""

from typing import Literal, TypeVar
from Vertex import Vertex

EdgeType = TypeVar('EdgeType', bound='Edge')

class Edge(object):

    def __init__(self, 
        vertices: "list[Vertex]",
        bond_type: Literal["", "=", "#"], 
        index: int
    ) -> None:
        self.index: int = index
        self.atoms: "list[Vertex]" = vertices
        self.bond_type: Literal["", "=", "#"] = bond_type
        self.symbols: list[str] = list(str(vertex.symbol) for vertex in vertices)
        self.indices: list[int] = list(int(vertex.index) for vertex in vertices)
        self.core_type: bool = not 'R' in self.symbols

    def __eq__(self: EdgeType, __o: EdgeType) -> bool:
        # Validates structural edge equality
        # self === the functional group
        # __o === target organic molecule

        # core_type edge matching
        if (self.core_type):
            return (self.bond_type == __o.bond_type) and (set(self.symbols) == set(__o.symbols))

        # R_type edge matching
        return (self.bond_type == __o.bond_type) and ([v for v  in self.symbols if v != 'R'][0] in __o.symbols)
    
    def __hash__(self):
        return hash((self.bond_type, self.symbols[0], self.symbols[1], self.core_type))

    def __str__(self):
        return self.bond_type.join([''.join([str(r) for r in list(z)]) for z in zip(self.symbols, self.indices)])

    def __repr__(self):
        return self.bond_type.join([''.join([str(r) for r in list(z)]) for z in zip(self.symbols, self.indices)])

    def mirror(self: EdgeType, __o: EdgeType):
        # Validates a mirrored edge equality (i.e if two edge objects are the exact same edge in a particular graph)
        return (self.bond_type == __o.bond_type) & (self.symbols == __o.symbols) & (self.indices == __o.indices)

    def complement_vertex(self, vertex_index):
        # gets the complement vertex in the edge when given a specific vertex_index
        return [vertex for vertex in self.atoms if vertex.index != vertex_index][0]
        # try:
        # except:
        #     ValueError("There was an error in the complement calculations")