"""The Vertex is representative of an atom in a molecule, and it is a vertex of part of the simple undirected connected molecular graph (graph theory)

The Vertex index is the particular index of an atom inside of the SMILES code
The Vertex symbol is the atomic symbol associated with a particular atom at a particular index in the SMILES code
The Vertex ring_type is the classification of the atom within a ring
(Atom class)
"""

from typing import Literal

class Vertex():

    def __init__(self, index: int = 0, symbol: str = 'C', isAromatic: bool = False, valenceElectronsRequired: int = 0, charge: str = ""):
        self.index: int = index
        self.symbol: str = symbol
        self.ring_type: Literal["aromatic", "non-aromatic", "non-cyclic"] = "aromatic" if isAromatic else "non-cyclic"
        self.implicit_degree: int = 0               # the number of single bonded hydrogens
        self.explicit_degree: int = 0               # the number of explicitly determined edges from the SMILES code
        self.total_degree: int = self.explicit_degree + self.implicit_degree # the sum of the explicit and implicit degrees
        self.valence_electrons_required: int = valenceElectronsRequired
        self.charge = charge

        ##### Charge Valence Electrons #####
        if charge == '+': self.valence_electrons_required+=1
        elif charge == '-': self.valence_electrons_required-=1

    def __eq__(self, value):
        """ Same symbol and index imply equivalent atoms """

        return self.index == value.index and self.symbol == value.symbol

    def __str__(self):
        return str(self.symbol) + str(self.index)

    def __repr__(self):
        return str(self.symbol) + str(self.index)
