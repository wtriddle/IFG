"""A vertex of a molecular graph, representative of an atom in a molecule"""

from typing import Literal

class Vertex():
    """A vertex of a molecular graph.

        Parameters
        ----------
        index : int
            The unique integer index identifier for a particular vertex in a molecular graph
        
        symbol : str
            The atomic symbol inclusive of charge for the vertex object

        is_aromatic : bool
            The assertion of if the vertex is aromatic or non-aromatic

        valence_electrons_required : int
            The integer number of valence electrons required to fill the outermost atomic orbital, considering impact of charge

        charge: str
            The symbol of a charge attached to a vertex

    
        Returns 
        -------
        Vertex
            The vertex object representative of an atom in a molecular graph
    
    """

    def __init__(self, 
        index: int = 0, 
        symbol: str = 'C', 
        is_aromatic: bool = False, 
        valence_electrons_required: int = 0, 
        charge: str = ""
    ):
        """Generates a new Vertex for a software molecule graph"""

        ##### Vertex Identifiers #####
        self.index: int = index
        """The unique integer index identifier for this vertex"""

        self.symbol: str = symbol
        """The atomic symbol inclusive of charge for this vertex"""

        ##### Ring Classification #####
        self.ring_type: Literal["aromatic", "non-aromatic", "non-cyclic"] = "aromatic" if is_aromatic else "non-cyclic"
        """The type of ring, if any, that a vertex is apart of"""

        ##### Degree & Atomic Valence Structure #####
        self.implicit_degree: int = 0
        """The number of hidden hydrogens connected to the vertex"""

        self.explicit_degree: int = 0
        """The number of explictly connected edges to the vertex"""

        self.total_degree: int = self.explicit_degree + self.implicit_degree
        """The total number of bonds (single, double, or triple) connected to a vertex"""

        self.valence_electrons_required: int = valence_electrons_required
        """The total number of valence electrons required to fill the outermost atomic orbital (i.e. its preferred electron configuration)"""
        
        self.charge: str = charge
        """The charge symbol attached to an atomic vertex"""

        ##### Charge Valence Electrons #####
        if charge == '+': self.valence_electrons_required+=1
        elif charge == '-': self.valence_electrons_required-=1

    def __str__(self):
        """Vertex String Representation"""
        return str(self.symbol) + str(self.index)

    def __repr__(self):
        """Vertex General Representation (Same as String Representation)"""
        return str(self.symbol) + str(self.index)
