"""The Vertex is representative of an atom in a molecule, and it is a vertex of part of the simple undirected connected molecular graph (graph theory)

The Vertex index is the particular index of an atom inside of the SMILES code
The Vertex symbol is the atomic symbol associated with a particular atom at a particular index in the SMILES code
The Vertex ring_type is the classification of the atom within a ring
The Vertex implicit_degree is the number of single bonded hydrogens implicitly connected to a vertex to satisfy its valence structure
The Vertex explicit_degree is the number of edges where the vertex appears in the edge set
The Vertex total_degree is the sum of the implicit and explicit degrees
The Vertex valence_electrons_required is the number of electrons required for a Vertex (atom) to fill its atomic orbital, charge inclusive
The Vertex charge is the assigned charge symbol of a Vertex and alters the valence_electrons_required 
"""

from typing import Literal

class Vertex():
    """A Vertex of a software molecule graph"""

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
        self.symbol: str = symbol

        ##### Ring Classification #####
        self.ring_type: Literal["aromatic", "non-aromatic", "non-cyclic"] = "aromatic" if is_aromatic else "non-cyclic"

        ##### Degree & Atomic Valence Structure #####
        self.implicit_degree: int = 0
        self.explicit_degree: int = 0
        self.total_degree: int = self.explicit_degree + self.implicit_degree
        self.valence_electrons_required: int = valence_electrons_required
        self.charge = charge

        ##### Charge Valence Electrons #####
        if charge == '+': self.valence_electrons_required+=1
        elif charge == '-': self.valence_electrons_required-=1

    def __str__(self):
        """Vertex String Representation"""
        return str(self.symbol) + str(self.index)

    def __repr__(self):
        """Vertex General Representation (Same as String Representation)"""
        return str(self.symbol) + str(self.index)
