""" Represents a SMILES code as its molecular equivalent with atomic properties.

A molecule class contains various containers, holding numerical data, string data, and Atom objects
It derives each all atoms and their bonds entirley from a string SMILEScode, which contains the
ingredients for how a molecule is connected and bonded in space. This Molecule class breaks down those
ingredients to extraplote basic bonding and molecular data which comprise the molecule described by the SMILES.

This molecule class only supports CHON sets.

Basic usage:

    // For a smiles code
    mol = Molecule('O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O', 'ABEGOH')
    
    // For a functioanl group template
    functionalGroup = Molecule('RC(=O)O', 'Carboxylicacid')

Key Attributes:
    atomData: Dictionary with atomIndex key to Atom Object value pairing
    bondData: Dictionary with atomIndex key to list of Atom Object value pairing
    ringData: Dictionary with ring key to count of specific ring type value

Notes:
    SMILES codes are exported from this module in all upper case
    The module contains Atom objects which are simply symbol index objects 
"""

from Atom import Atom
import re
from collections import OrderedDict


class Molecule():
    """ Object representation of a SMILES code based on symbol by symbol analysis 

    Basic usage:

    >> mol = Molecule('O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O', 'ABEGOH')
    >> print(mol)
    ABEGOH : O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O
    >> functionalGroup = Molecule('RC(=O)O', 'Carboxylicacid')
    """

    def __init__(self, smiles, name):
        """ Intializes molecule obejct with atom and bonding data for all symbols
            Derives Cyclic and noncyclic properties for all atoms
            Derives ring counts
            Contains symbol matching lists for string chars
            Derives alcohol counts, inclucing cyclic/aromatic attacments

        """
        # Input data
        self.SMILES = self.formatSmiles(smiles)
        self.NAME = name

        # Symbol matching lists and regex's
        self.ATOM_REGEX = re.compile(r'[a-zA-Z]')
        self.CHARGE_REGEX = re.compile(r'\+|\-')
        self.BOND_REGEX = re.compile(r'\=|\#')
        self.ATOMS = ['C', 'O', 'N', 'X', 'Z', 'S', 'I',
                      'F', 'c', 'n', 'o', 'x', 'z', 's', 'i', 'f', 'R', 'P', 'p']
        self.BONDS = ['=', '#']
        self.BRACKETS = ['[', ']']
        self.CHARGES = ['+', '-']
        self.NUMBERS = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
        self.LINEARSYMBOLS = self.ATOMS + self.BONDS + self.BRACKETS + self.CHARGES

        # Ring Containers
        self.RING_SELF_NUM_POSITIONS = []
        self.RING_SELF = {}
        self.RING_PARTNER_NUM_POSITIONS = []
        self.RING_PARTNERS = {}

        self._RING()  # Initalize Ring Containers

        # Index Containers
        self.AROMATICINDICES = []
        self.CYCLICINDICES = []

        self._INCIDES()  # Initalize Index Containers

        # Collect ring data from ring containers
        self.ringCount = len(self.RING_SELF) / 2
        self.aromaticCount = 0
        self.nonAromaticCount = 0
        self._RING_COUNTS()

        # Dictionary mapping counts to keys
        self.ringData = {
            "aromaticRingCount": self.aromaticCount,
            "nonAromaticRingCount": self.nonAromaticCount,
            "ringCount": self.ringCount}

        # Once aromaticity is extracted, analysis of SMILES codes can be done in upper case letters only for simplicity
        self.SMILES = self.SMILES.upper()
        self.ALCOHOLICINDICES = []
        self.atomData = self.initializeAtomData()
        self.atomCount = len(self.atomData)
        self.bondData = self.initializeBondData()
        self.chargedMol = True if len(
            self.CHARGE_REGEX.findall(smiles)) != 0 else False
        self.AMINOACID = True if len(re.compile(
            r'\[[nN]H[23]?\+\]').findall(smiles)) != 0 else False

    def __str__(self):
        return self.NAME + " : " + self.SMILES

    def getSymbolDict(self):
        """ Create a dictionary of index to atom symbol of the molecule's atoms based on its current state"""
        symbolDict = {
            atom.index: atom.symbol
            for atom in self.atomData.values()
        }
        return symbolDict

    def initializeAtomData(self):
        """ Return dictionary of atomData based on all symbols in SMILES code.

            Atom objects are created based on the symbol and index position of that atom
            Handles charges and finds valid alcohol indices
        """

        atomData = {}
        atomIndex = -1

        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1
                atomSymbol = symbol

                # If atom charged, collect its bracketed group
                if self.SMILES[pos-1] == '[':
                    atomSymbol = self.getChargedGroup(pos-1)

                atom = Atom(atomIndex, atomSymbol)
                atomData[atomIndex] = atom

                # Alochols stem from oxygens
                if symbol == 'O':
                    self.determineAlcoholGroup(pos, atomIndex)

        return atomData

    def initializeBondData(self):
        """ Returns the bondData dictionary after SMILES code analysis.

            From a given atom position, the left and right bond datas are collected with their respective algorithms
            Each bond represents a single Atom object. However, bondData keys to a list of Atom objects, representing 
            all atoms bonded to the atom of interest.

        """

        bondData = {}
        atomIndex = -1

        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1
                bonds = []

                # Deal with charges based on special starting positions
                if self.SMILES[pos-1] == '[':
                    leftBond = self.getLeftBond(LNPos=pos-1, LNIndex=atomIndex)
                    rightBonds = self.getRightBonds(
                        RNPos=pos+2, RNIndex=atomIndex)
                else:
                    leftBond = self.getLeftBond(LNPos=pos, LNIndex=atomIndex)
                    rightBonds = self.getRightBonds(
                        RNPos=pos, RNIndex=atomIndex)

                if leftBond:
                    bonds.append(leftBond)

                if rightBonds:
                    for rightBond in rightBonds:
                        bonds.append(rightBond)

                bondData[atomIndex] = bonds
                del(bonds)

        return bondData

    def getLeftBond(self, LNPos, LNIndex):
        """ Returns an single Atom object. 
            Return an empty array [] if no LN exists

            LNPos (int) : string position of the atom whose left bond is to be found
            LNIndex (int) : index position of that atom within the smiles code

            Notes:
                LN stands for Left Neighbor
        """

        # Begin one position to left
        LNPos -= 1
        if LNPos < 0:
            return []
        LN = self.SMILES[LNPos]

        # An explicit bond exists for a left-bonded atom if and only if the bond is located one position to the left of the atom
        explicitBond = LN if LN in self.BONDS else ""
        if explicitBond:
            LNPos -= 1
            LN = self.SMILES[LNPos]

        scope = 0  # Parenthesis depth of nesting value. Initally 0, meaning to look for neighbors on its own layer or nested scope
        leftBond = []

        # Loop continues until a LINEARSYMBOL is found on the same scope as its bonded neighbor
        while LN not in self.LINEARSYMBOLS or scope > 0:

            # Special case to pass double parenthesis left neighbors, ie. N(C)(C).
            # The second parenthesied C needs to find the N as its proper LN
            # This statement allows passage through its other C neighbor, which it is NOT directly attachted to
            if LN == ')' and scope == -1:
                scope = 0

            if LN == ')':
                scope += 1
            elif LN == '(':
                scope -= 1
            elif LN in self.ATOMS:
                LNIndex -= 1

            LNPos -= 1
            if LNPos < 0:
                return []
            LN = self.SMILES[LNPos]

        else:
            LNIndex -= 1

            # If LN is a charged group, overwrite it
            if LN == ']':
                LN = self.getChargedGroup(LNPos)

            leftBond = Atom(LNIndex, explicitBond + LN)

        return leftBond

    def getRightBonds(self, RNPos, RNIndex):
        """ Returns an atom list of all right bonds for a given atom.
            Returns an empty list if no right bond exists


            RNPos (int) : The position of the atom in the SMILES code whose right bonds are to be found
            RNIndex (int) : The index position of the atom in the SMILES code

            Notes:
                This is a recursive function
                RN stands for Right Neighbor
                Right groups split to different potions of the smiles code (dynamic symbols)
                The two dynamic symbols are () and numbers. Both are handled accordingly here. 
        """

        # Begin at one position to right
        RNPos += 1
        if RNPos >= len(self.SMILES):
            return []
        RN = self.SMILES[RNPos]

        # A closing parenthesis explicitly means no right bonds
        if RN == ')':
            return []

        scope = 0  # Parenthesis depth of nesting value. Initally 0, meaning to look for neighbors on its own layer or nested scope
        rightBonds = []

        while RN not in self.LINEARSYMBOLS or scope > 0:

            if RN == '(':

                # Recursion to find first nested bond. Must be on same scope as atom (i.e scope 0)
                innerParenthBond = self.getRightBonds(
                    RNPos, RNIndex
                ) if scope == 0 else ''

                # Function returns list, capture first and only atom element
                if innerParenthBond:
                    rightBonds.append(innerParenthBond[0])

                # Increment scope, now going inside a parenthesis
                scope += 1

            if RN == ')':
                scope -= 1
            if RN in self.ATOMS:
                RNIndex += 1
            if scope == 0 and RN in self.NUMBERS:
                # Every number always has a bonded partner
                numGroup = self.numbersHandler(RNPos)
                rightBonds.append(numGroup)

            RNPos += 1
            if RNPos >= len(self.SMILES) or scope < 0:
                return rightBonds
            RN = self.SMILES[RNPos]

        else:
            RNIndex += 1

            if RN in self.ATOMS:
                rightBonds.append(Atom(RNIndex, RN))

            elif RN in self.BONDS:
                explicitBond = RN
                RNPos += 1
                RN = self.SMILES[RNPos]

                if self.SMILES[RNPos] == '[':
                    RN = self.getChargedGroup(RNPos)

                rightBonds.append(Atom(RNIndex, explicitBond + RN))

            elif RN == '[':
                chargedGroup = self.getChargedGroup(RNPos)
                rightBonds.append(Atom(RNIndex, chargedGroup))

        return rightBonds

    def getChargedGroup(self, pos):
        """ Returns a charged group based on a selected bracket position in the its smiles code
            pos (int): string position of an openeing or closing bracket in the smiles code
        """

        # String that holds all characters that make up the bracketed charge group
        chargedGroup = ''

        # Reverse case
        if self.SMILES[pos] == ']':
            while self.SMILES[pos] != '[':
                chargedGroup += self.SMILES[pos]
                pos -= 1
            return '[' + chargedGroup[::-1]

        # Forwards case
        if self.SMILES[pos] == '[':
            while self.SMILES[pos] != ']':
                chargedGroup += self.SMILES[pos]
                pos += 1
            return chargedGroup + ']'

        # Capture errors in positions
        if self.SMILES[pos] not in self.BRACKETS:
            raise ValueError(
                "The position passed does not point to a bracket in the smiles code")

    def _RING(self):
        """ Initalizes the four ring data containers
            RING_SELF: number position to atom data dict
            RING_PARTNERS: number position to complementary number data
            RING_SELF_NUM_POSITIONS: positions of the numbers which began the partner relation
            RING_PARTNER_NUM_POSITIONS: positions of the numbers which ended the partner relation


            Notes: 
                evalutedNumbers is a dictionary of number symbol : position in smiles code
        """
        atomIndex = -1
        evaluatedNumbers = {}
        atomSymbol = ''

        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1
                atomSymbol = symbol

            if symbol in self.NUMBERS:

                if self.SMILES[pos-1] == ']':
                    atomSymbol = self.getChargedGroup(pos-1)

                atom = Atom(atomIndex, atomSymbol)

                # Holds data for num position to imediate bond neighbor pairing, i.e. C1, =N2, O3, etc.
                self.RING_SELF[pos] = atom

                # If the partner number has been found, then the bonded atom has been found
                if symbol in evaluatedNumbers:

                    # Get inital number position, where the ring began
                    initalNumberPos = evaluatedNumbers[symbol]

                    # Update the partner dict to hold position : connected atom info
                    self.RING_PARTNERS[initalNumberPos] = atom
                    self.RING_PARTNERS[pos] = self.RING_SELF[initalNumberPos]

                    # Multiple of the same number can occur in a SMILES code, so open that number path for another pairing
                    del(evaluatedNumbers[symbol])

                    # Position of where the rings end in the SMILES code
                    self.RING_PARTNER_NUM_POSITIONS.append(pos)
                    continue

                # Position of where the rings start in the SMILES code
                self.RING_SELF_NUM_POSITIONS.append(pos)

                # Number to positioning key:value pairing
                evaluatedNumbers[symbol] = pos

    def numbersHandler(self, pos):
        """ Return the partner of a given position to retrieve what atom a specific number is connected to

            pos (int) : Position of a ring whose partner atom is to be determined

        """
        return self.RING_PARTNERS[pos]

    def determineAlcoholGroup(self, pos, atomIndex):
        """ Updates the Alochol Index list with the index of the oxygen involved in a valid alcohol 

            pos (int) : position of the oxygen to be checked for alocholic properties
            atomIndex (int) : index position of said oxygen

            Notes:
                Only four cases where an alcohol is determined. 
                Currently alochol cases are hard coded here for CHON set of data.
        """

        # Final atom alcohol group
        if pos == len(self.SMILES) - 1:
            if self.SMILES[pos-1].upper() == 'C' or self.SMILES[pos-1] in self.NUMBERS or self.SMILES[pos-1] == ')':
                self.ALCOHOLICINDICES.append(atomIndex)

        # First atom alcohol group
        elif atomIndex == 0 and self.SMILES[pos+1].upper() == 'C':
            self.ALCOHOLICINDICES.append(atomIndex)

        # Lone alcohol group
        elif self.SMILES[pos-1:pos+2] == '(O)':
            self.ALCOHOLICINDICES.append(atomIndex)

        # Ending parenthesis alcohol group
        elif self.SMILES[pos+1] == ')' and self.SMILES[pos-1] not in self.BONDS and self.SMILES[pos-1] not in self.BRACKETS:
            self.ALCOHOLICINDICES.append(atomIndex)

        return 0

    def _RING_COUNTS(self):
        """ Loops over Ring containers to determine the counts of aromatic and non aromatic rings in the molecule

            Notes:
                Uses the two atoms involved in the number ring positioning to detemine what kind of ring the atoms are apart of
                Some extra neighbor checking is required to confirm aromatic/non aromatic presence fully
        """

        # Simplification of aromatic/nonaromatic count
        allAtoms = ''.join(self.ATOM_REGEX.findall(self.SMILES))

        if allAtoms.islower():
            self.aromaticCount = self.ringCount
            return 0

        elif allAtoms.isupper():
            self.nonAromaticCount = self.ringCount
            return 0

        # Upper is nonaromatic, lower is aromatic
        for pos in self.RING_SELF_NUM_POSITIONS:
            ring = self.RING_SELF[pos]
            ringPartner = self.RING_PARTNERS[pos]

            # Both are aromatic
            if ring.symbol.islower() and ringPartner.symbol.islower():

                if self.SMILES[pos+1] in self.ATOMS:

                    if self.SMILES[pos+1].islower():
                        self.aromaticCount += 1

                    elif self.SMILES[pos+1].isupper():
                        self.nonAromaticCount += 1

                else:
                    self.aromaticCount += 1

            # Both are non aromatic
            elif ring.symbol.isupper() and ringPartner.symbol.isupper():

                if self.SMILES[pos+1] in self.ATOMS:

                    if self.SMILES[pos+1].islower():
                        self.aromaticCount += 1

                    elif self.SMILES[pos+1].isupper():
                        self.nonAromaticCount += 1

                else:
                    self.nonAromaticCount += 1

            # They are differnt, must be nonaromatic
            else:
                self.nonAromaticCount += 1

    def _INCIDES(self):
        """ Determines the cyclic/aromatic atom indicies from the smiles code.

            The indices are placed into pre-intialized self containers, see the __init__ call order

            Notes:
                Scope and depth are used. 
                If a ring ceases inside a parenthesis group, then all indicies in that scope and above it are contained in the ring
                Therefore, tracking of indices on specific scope levels (parenthesis scopes) is necessary to obtain the most accurate
                ring indices
        """

        evaluatedNumbers = []

        atomIndex = -1
        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1

            if symbol.islower() and atomIndex not in self.AROMATICINDICES:
                self.AROMATICINDICES.append(atomIndex)
                continue

            # If a new ring has been found in the SMILES code via a number
            if symbol in self.NUMBERS and pos not in evaluatedNumbers:

                evaluatedNumbers.append(pos)

                # Index by how deeply nested a parenthesis group is.
                # I.e. 0th depth is scope where number was found,
                # 1st is an inner parethensis (), 2nd is (... (...) ... ), 3rd is (... (... (...) ...) ...), etc.
                scopeIndices = [[]]
                scope = 0
                scopeIndices[0].append(atomIndex)
                RNPos = pos + 1
                RNindex = atomIndex
                RN = self.SMILES[RNPos]

                # Symbol is equal to SMILES number which opened the ring. Loop untl the same number is encountered by RN
                while RN != symbol:

                    if RN in self.ATOMS:
                        RNindex += 1
                        # If a new nested parenthesis or deeper parenthesis is found, i.e. in between a parenthesis group (...()...)
                        if len(scopeIndices) == scope:
                            # Then create a new list tracking the inner indices of that nested parenthesis group
                            scopeIndices.append([RNindex])
                        else:  # If the parenthesis depth did not go any deeper, then access the current depth and append the atom index to it
                            scopeIndices[scope].append(RNindex)

                    if RN == '(':
                        scope += 1
                    if RN == ')':
                        # A closing parenthesis means the scope has risen, and that the number cannot exist within that specifc nested depth scope
                        # Remove the tracked indices, as they no longer can be cyclic relative to the current scope that the numbered ring being tracked is on
                        del(scopeIndices[scope])
                        scope -= 1
                    RNPos += 1
                    RN = self.SMILES[RNPos]

                else:
                    evaluatedNumbers.append(RNPos)
                    # If a number ends within a depth deeper than 0, i.e. ...1... (... (...1)...),
                    # Then all all indices in in between must be cyclic
                    # Some nested scopes may go deeper, but may not conclude the ring.
                    # For example, ...1... (... (...) ... 1). The middle parenthesis is not apart of the ring structure related to the number 1.
                    # However, the number 1 still closes within a nested parenthesis. Because the depth of conclusion is not known before this algorithm is run
                    # All depths within indivudal depth counts must be tracked to obtain accurate cyclic index information
                    for scope in scopeIndices:
                        for index in scope:  # Add all indices within that scope to self.CYCLICINDICES
                            if index not in self.CYCLICINDICES:
                                self.CYCLICINDICES.append(index)

    def formatSmiles(self, smiles):
        """ Remove [H+] symbols entirley from a smiles code and return a DLA-SAR converted smiles code """
        reFormatted = ""
        for pos, symbol in enumerate(smiles):

            if symbol == '[':
                startBracketPos = pos

            if symbol == 'H':
                cutPos = pos
                while smiles[cutPos] != ']':
                    cutPos += 1
                reFormatted = smiles[0:startBracketPos] + \
                    smiles[startBracketPos+1] + smiles[cutPos+1:len(smiles)]
                reFormatted = self.formatSmiles(reFormatted)
                break

        if reFormatted == "":
            return self.DLAtoSARconversion(smiles)
        else:
            return self.DLAtoSARconversion(reFormatted)

    def DLAtoSARconversion(self, template):
        """ Tranform double letterd atoms (DLA) to single atom representations (SAR)
            Required to perform sybmol by symbol analysis on the smiles code

            Notes:
                Placeholder for futher analysis upon non CHON atom. This is untested at large, program can only handle CHONS at the moment
        """
        # Conversion Legend for DLA's to SAR's
        legend = {
            "Br": "X",
            "Cl": "Z"
        }

        pos = -1
        newTemplate = ""

        while pos != len(template) - 1:

            pos += 1
            # Potential DLA from two letters in SMILES code
            DLA = template[pos:pos+2]

            # Convert DLA'S via the legend
            if DLA in legend:
                newTemplate += legend[DLA]
                pos += 1

            # If no DLA is found, keep the template the same
            else:
                newTemplate += template[pos]

        return newTemplate
