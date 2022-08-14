"""An Edge is representative of a bond in the molecule, and appears as the edges of a simple undirected connected graph (graph theory)

Edge bond_type is the type of bond between two verticies (i.e atoms) in the molecular graph
Edge verticies are the two atom symbols which are connected by the edge
Edge indices are the two atom index values which identify the two atoms that are connected by the edge (local to the graph)
Edge core_type is the assertion of both verticies being non-R atoms (the negation of core_type is called an R_type which has an R atom involved in the edge connection)
"""

from typing import Literal, TypeVar
from Vertex import Vertex

EdgeType = TypeVar('EdgeType', bound='Edge')

class Edge(object):

    def __init__(self, 
        atoms: "list[Vertex]",
        bond_type: Literal["", "=", "#"], 
    ) -> None:
        self.bond_type: Literal["", "=", "#"] = bond_type
        self.vertices: list[str] = list(str(atom.symbol) for atom in atoms)
        self.indices: list[int] = list(int(atom.index) for atom in atoms)
        self.core_type: bool = not 'R' in self.vertices

    def __eq__(self: EdgeType, __o: EdgeType) -> bool:
        # Validates structural edge equality
        # self === the functional group
        # __o === target organic molecule

        # core_type edge matching
        if (self.core_type):
            return (self.bond_type == __o.bond_type) & (self.vertices == __o.vertices)

        # R_type edge matching
        return (self.bond_type == __o.bond_type) & (self.vertices.remove('R') in __o.vertices)

    def __str__(self):
        return self.bond_type.join([''.join([str(r) for r in list(z)]) for z in zip(self.vertices, self.indices)])

    def __repr__(self):
        return self.bond_type.join([''.join([str(r) for r in list(z)]) for z in zip(self.vertices, self.indices)])

    def mirror(self: EdgeType, __o: EdgeType):
        # Validates a mirrored edge equality (i.e if two edge objects are the exact same edge in a particular graph)
        return (self.bond_type == __o.bond_type) & (self.vertices == __o.vertices) & (self.indices == __o.indices)
