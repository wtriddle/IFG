"""The Vertex is representative of an atom in a molecule, and it is a vertex of part of the simple undirected connected molecular graph (graph theory)

The Vertex index is the particular index of an atom inside of the SMILES code
The Vertex symbol is the atomic symbol associated with a particular atom at a particular index in the SMILES code
The Vertex ring_type is the classification of the atom within a ring
(Atom class)
"""

from typing import Literal

class Vertex():

    def __init__(self, index: int, symbol: str, isAromatic: bool = False):
        self.index: int = index
        self.symbol: str = symbol
        self.ring_type: Literal["aromatic", "non-aromatic", "non-cyclic"] = "aromatic" if isAromatic else "non-cyclic"

    def __eq__(self, value):
        """ Same symbol and index imply equivalent atoms """

        return self.index == value.index and self.symbol == value.symbol

    def __str__(self):
        return str(self.symbol) + str(self.index)

    def __repr__(self):
        return str(self.symbol) + str(self.index)
