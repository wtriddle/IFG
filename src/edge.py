"""An Edge is representative of a bond in the molecule, and appears as the edges of a simple undirected connected graph (graph theory)

The Edge index is the index of the edge in the graph
The Edge atoms are the Vertex objects involved in the edge
The Edge indices are the two vertex indices involved in the edge (these indices are local to the graph and cannot be overwritten)
The Edge bond_type is the type of bond between two vertices (i.e atoms) in the molecular graph
The Edge symbols are the two atomic symbols of the vertices involved in the edge
The Edge core_type is if both symbols of the edge are non-R atoms (the negation of core_type is called an R_type which has an R atom involved in the edge connection)
"""

from typing import Literal, TypeVar
from vertex import Vertex

##### Edge Class Type #####
EdgeType = TypeVar('EdgeType', bound='Edge')

class Edge():
    """An Edge that connects two vertices of a software molecule graph together"""

    def __init__(self, 
        vertices: "list[Vertex]",
        bond_type: Literal["", "=", "#"], 
        index: int
    ) -> None:
        """Generates a new Edge for a software molecule graph"""

        ##### Unique Edge Identifier #####
        self.index: int = index

        ##### Vertex Identifiers #####
        self.atoms: "list[Vertex]" = vertices
        self.indices: list[int] = list(int(vertex.index) for vertex in vertices)

        ##### Structural Edge Identifiers #####
        self.bond_type: Literal["", "=", "#"] = bond_type
        self.symbols: list[str] = list(str(vertex.symbol) for vertex in vertices)
        self.core_type: bool = not 'R' in self.symbols
        

    def __eq__(self: EdgeType, __o: EdgeType) -> bool:
        """Tests Structural Edge Equality"""
        return (self.bond_type == __o.bond_type) and (set(self.symbols) == set(__o.symbols))
    
    def __hash__(self):
        """Hash of the Edge defined by the Structural Identifier"""
        return hash((self.bond_type, self.symbols[0], self.symbols[1], self.core_type))

    def __str__(self) -> str:
        """String Representation of the Edge"""
        return self.bond_type.join([''.join([str(r) for r in list(z)]) for z in zip(self.symbols, self.indices)])

    def __repr__(self) -> str:
        """General Representation of the Edge (Equivalent to its String Representation)"""
        return self.bond_type.join([''.join([str(r) for r in list(z)]) for z in zip(self.symbols, self.indices)])

    def mirror(self: EdgeType, __o: EdgeType) -> bool:
        """Tests Unique Edge (Inclusive of Strucutre) Equality"""
        return (self.bond_type == __o.bond_type) & (self.symbols == __o.symbols) & (self.index == __o.index)

    def complement_vertex(self, vertex_index) -> Vertex:
        """Gets the complementary (or opposite) vertex involved in the edge when given one of the two vertex indices of the edge"""
        return [vertex for vertex in self.atoms if vertex.index != vertex_index][0]