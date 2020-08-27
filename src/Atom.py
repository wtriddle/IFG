class Atom():

    def __init__(self, index, atomSymbol):
        self.index = index
        self.symbol = atomSymbol

    def __eq__(self, value):
        return self.index == value.index and self.symbol == value.symbol

    def __str__(self):
        return str(self.symbol)

    def __repr__(self):
        return self.symbol
