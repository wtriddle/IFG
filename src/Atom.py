""" Represents a single atom part of any larger structure.

Key Attributes:
    Index: The index location of a specific atom within a larger molecule structure, such as the SMILES code of functional group
    Symbol: The string representing the atom. It could be a charged group and be more than just a single symbol
"""


class Atom():
    """ Represents a symbol and index as a self-contained Atom object """

    def __init__(self, index, atomSymbol):
        self.index = index
        self.symbol = atomSymbol

    def __eq__(self, value):
        """ Same symbol and index imply equivalent atoms """

        return self.index == value.index and self.symbol == value.symbol

    def __str__(self):
        return str(self.symbol)

    def __repr__(self):
        return self.symbol
