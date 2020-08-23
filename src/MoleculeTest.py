''' This file will be altered and checked against Molecule after logic changes '''

import re


class MoleculeTest():
    ''' Object representation of a SMILES code based on symbol by symbol analysis'''

    def __init__(self, smiles, name):

        # Input data
        self.SMILES = self.formatSmiles(smiles)
        self.NAME = name

        # Symbol matching lists and regex's
        self.ATOM_REGEX = re.compile(r'[a-zA-Z]')
        self.CHARGE_REGEX = re.compile(r'\+|\-')
        self.BOND_REGEX = re.compile(r'\=|\#')
        self.ATOMS = ['C', 'O', 'N', 'X', 'Z', 'S', 'I',
                      'F', 'c', 'n', 'o', 'x', 'z', 's', 'i', 'f', 'R']
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
        self.RINGDICT = {
            "aromaticRings": self.aromaticCount,
            "nonAromaticRings": self.nonAromaticCount,
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

    def initializeAtomData(self):

        atomData = []
        atomIndex = -1

        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1

                if self.SMILES[pos-1] == '[':
                    chargedGroup = self.SMILES[pos-1:pos+3]
                    atom = [atomIndex, chargedGroup]
                else:
                    atom = [atomIndex, symbol]
                atomData.append(atom)

                if symbol == 'O':
                    self.determineAlcoholGroup(pos, atomIndex)

        return atomData

    def initializeBondData(self):

        bondData = []
        atomIndex = -1

        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1
                bonds = []

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

                bondData.append(bonds)
                del(bonds)

        return bondData

    def getLeftBond(self, LNPos, LNIndex):

        LNPos -= 1
        if LNPos < 0:
            return 0
        LN = self.SMILES[LNPos]

        expBond = ""
        if LN in self.BONDS:
            expBond = LN
            LNPos -= 1
            LN = self.SMILES[LNPos]

        scope = 0
        leftBond = []

        while LN not in self.LINEARSYMBOLS or scope > 0:

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
                return 0
            LN = self.SMILES[LNPos]

        else:
            LNIndex -= 1

            if LN == ']':
                chargedGroup = ']'

                while self.SMILES[LNPos] != '[':
                    LNPos -= 1
                    chargedGroup = self.SMILES[LNPos] + chargedGroup

                if expBond:
                    leftBond = [LNIndex, expBond + chargedGroup]
                else:
                    leftBond = [LNIndex, chargedGroup]

            elif LN in self.ATOMS:

                if expBond:
                    leftBond = [LNIndex, expBond + LN]
                else:
                    leftBond = [LNIndex, LN]

        return leftBond

    def getRightBonds(self, RNPos, RNIndex):

        RNPos += 1
        if RNPos >= len(self.SMILES):
            return 0
        RN = self.SMILES[RNPos]

        if RN == ')':
            return 0

        scope = 0
        parenthGroup = ""
        parenthGroups = []
        rightBonds = []

        while RN not in self.LINEARSYMBOLS or scope > 0:

            if RN == '(':
                scope += 1
            if scope > 0 and not parenthGroup:
                parenthGroups.append([RNIndex+1])
            if scope > 0:
                parenthGroup += RN
            if RN == ')':
                scope -= 1
            if scope == 0 and parenthGroup:
                parenthGroups[-1].append(parenthGroup)
            if RN in self.ATOMS:
                RNIndex += 1
            if scope == 0 and RN in self.NUMBERS:
                numGroup = self.numbersHandler(RNPos)
                rightBonds.append(numGroup)
            RNPos += 1
            if RNPos >= len(self.SMILES) or scope < 0:
                break
            RN = self.SMILES[RNPos]

        else:
            RNIndex += 1

            if RN in self.ATOMS:
                rightBonds.append([RNIndex, RN])

            elif RN in self.BONDS:
                expBond = RN
                RNPos += 1

                if self.SMILES[RNPos] == '[':
                    chargedGroup = "["

                    while self.SMILES[RNPos] != ']':
                        RNPos += 1
                        chargedGroup += self.SMILES[RNPos]
                    rightBonds.append([RNIndex, expBond + chargedGroup])

                elif self.SMILES[RNPos] in self.ATOMS:
                    rightBonds.append([RNIndex, expBond + self.SMILES[RNPos]])

            elif RN == '[':
                chargedGroup = "["

                while self.SMILES[RNPos] != ']':
                    RNPos += 1
                    chargedGroup += self.SMILES[RNPos]
                rightBonds.append([RNIndex, chargedGroup])

            # parenthGroups in the form [first atom index in SMILESparenthesisGroup, SMILESparenthesisGroup]
            # SMILESparenthesisGroup includes both open and close parenthesis on both sides
            if parenthGroups:
                for group in parenthGroups:
                    pos = 1
                    symbol = group[1][pos]

                    if symbol in self.ATOMS:
                        rightBonds.append([group[0], symbol])

                    elif symbol in self.BONDS:
                        expBond = symbol
                        pos += 1

                        if group[1][pos] == '[':
                            chargedGroup = '['

                            while group[1][pos] != ']':
                                pos += 1
                                chargedGroup += group[1][pos]
                            rightBonds.append(
                                [group[0], expBond + chargedGroup])

                        elif group[1][pos] in self.ATOMS:
                            rightBonds.append(
                                [group[0], expBond + group[1][pos]])

                    elif symbol == '[':
                        chargedGroup = '['

                        while group[1][pos] != ']':
                            pos += 1
                            chargedGroup += group[1][pos]
                        rightBonds.append([group[0], chargedGroup])

        return rightBonds

    def getChargedGroup(self, pos, reverse=False):
        """ 

        Captures a charged group in reverse or forwards 

        """
        if reverse:
            chargedGroup = ''
            while self.SMILES[pos] != '[':
                pos -= 1
                chargedGroup += self.SMILES[pos]
            return chargedGroup[::-1]
        else:
            chargedGroup = '['
            while self.SMILES[pos] != '[':
                pos += 1
                chargedGroup += self.SMILES[pos]
            return chargedGroup

    def _RING(self):
        """ Initalizes the four ring data containers
            RING_SELF: number position to atom data dict
            RING_PARTNERS: number position to complementary number data
            RING_SELF_NUM_POSITIONS: positions of the numbers which began the partner relation
            RING_PARTNER_NUM_POSITIONS: positions of the numbers which ended the partner relation

        """
        atomIndex = -1
        evaluatedNumbers = {}
        atom = ''
        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1
                atom = symbol

            if symbol in self.NUMBERS:

                if self.SMILES[pos-1] == ']':
                    atom = self.getChargedGroup(pos, reverse=True)

                info = [atomIndex, atom]

                self.RING_SELF[pos] = info

                if symbol in evaluatedNumbers:
                    partnerPos = evaluatedNumbers[symbol]
                    self.RING_PARTNERS[partnerPos] = info
                    self.RING_PARTNERS[pos] = self.RING_SELF[partnerPos]
                    del(evaluatedNumbers[symbol])
                    self.RING_PARTNER_NUM_POSITIONS.append(pos)
                    continue

                self.RING_SELF_NUM_POSITIONS.append(pos)
                evaluatedNumbers[symbol] = pos

    def numbersHandler(self, pos):
        return self.RING_PARTNERS[pos]

    def determineAlcoholGroup(self, pos, atomIndex):

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
            if ring[1].islower() and ringPartner[1].islower():

                if self.SMILES[pos+1] in self.ATOMS:

                    if self.SMILES[pos+1].islower():
                        self.aromaticCount += 1

                    elif self.SMILES[pos+1].isupper():
                        self.nonAromaticCount += 1

                else:
                    self.aromaticCount += 1

            # Both are non aromatic
            elif ring[1].isupper() and ringPartner[1].isupper():

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
        """ Initalizes ring index lists for CYCLIC/AROMATIC indexed atoms """

        evaluatedNumbers = []
        self.CYCLICINDICES = []

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
                depthIndices = [[]]
                nestDepthCount = 0
                depthIndices[0].append(atomIndex)
                RNPos = pos + 1
                RNindex = atomIndex
                RN = self.SMILES[RNPos]

                # Symbol is equal to SMILES number which opened the ring. Loop untl the same number is encountered by RN
                while RN != symbol:

                    if RN in self.ATOMS:
                        RNindex += 1
                        # If a new nested parenthesis or deeper parenthesis is found, i.e. in between a parenthesis group (...()...)
                        if len(depthIndices) == nestDepthCount:
                            # Then create a new list tracking the inner indices of that nested parenthesis group
                            depthIndices.append([RNindex])
                        else:  # If the parenthesis depth did not go any deeper, then access the current depth and append the atom index to it
                            depthIndices[nestDepthCount].append(RNindex)

                    if RN == '(':
                        nestDepthCount += 1
                    if RN == ')':
                        # A closing parenthesis means the depth count has risen, and that the number cannot exist within that specifc nested depth scope
                        # Remove the tracked indies, as they no longer can be cyclic relative to the current numbered ring being tracked
                        del(depthIndices[nestDepthCount])
                        nestDepthCount -= 1
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
                    for scope in depthIndices:
                        for index in scope:  # Add all indices within that scope to self.CYCLICINDICES
                            if index not in self.CYCLICINDICES:
                                self.CYCLICINDICES.append(index)

    def formatSmiles(self, smiles):
        """ Remove [H+] symbols entirley from a smiles code and DLA-SAR convert the resulting string

                """
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
        ''' Tranform double letterd atoms (DLA) to single atom representations (SAR)
                        Required to perform sybmol by symbol analysis on the smiles code

                '''
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
