"""An edge of a molecular graph, representative of a bond in a molecule"""

from typing import Literal, TypeVar
from vertex import Vertex

##### Edge Class Type #####
EdgeType = TypeVar('EdgeType', bound='Edge')

class Edge():
    """An edge connecting two vertices in a molecular graph together.
    
        Parameters
        ----------
        vertices : "list[Vertex]"
            The two vertex objects connected together by an edge
        
        bond_type : Literal["", "=", "#"]
            The symbol of the bond involved in the molecular edge

        index: int
            The unique integer index identifier for a particular edge in a molecular graph
    
        Returns 
        -------
        Edge
            The edge object representative of a bond in a molecular graph
    """

    def __init__(self, 
        vertices: "list[Vertex]",
        bond_type: Literal["", "=", "#"], 
        index: int
    ) -> None:
        """Generates a new Edge for a software molecule graph"""

        ##### Unique Edge Identifier #####
        self.index: int = index
        """The index identifier of the edge"""

        ##### Vertex Identifiers #####
        self.atoms: "list[Vertex]" = vertices
        """The two vertex objects involved in the edge connection"""
        self.indices: list[int] = list(int(vertex.index) for vertex in vertices)
        """The two vertex indices involved in the edge connection"""

        ##### Structural Edge Identifiers #####
        self.bond_type: Literal["", "=", "#"] = bond_type
        """The type of bond between the two atomic vertices"""

        self.symbols: list[str] = list(str(vertex.symbol) for vertex in vertices)
        """The two symbols of the vertex objects involved in the edge connection"""

        self.core_type: bool = not 'R' in self.symbols
        """The assertion of a core type edge without any R vertices involved in the edge connection"""
        

    def __eq__(self: EdgeType, __o: EdgeType) -> bool:
        """Tests Structural Edge Equality.
        
            Structural equality is when the symbols of both vertices and the bond between them are both identical.

        """
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
        """Returns true if the index, bond and symbols involved in two different edges are all equivalent."""
        return (self.bond_type == __o.bond_type) & (self.symbols == __o.symbols) & (self.index == __o.index)

    def complement_vertex(self, vertex_index: int) -> Vertex:
        """Gets the complementary (or opposite) vertex involved in the edge when given one of the two vertex indices of the edge.

        Parameters
        ----------
        vertex_index: int
            The index of one of the vertices involved in the edge

        Returns
        -------
        Vertex
            The vertex object of the other vertex involved in the edge
        """
        return [vertex for vertex in self.atoms if vertex.index != vertex_index][0]