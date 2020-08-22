''' This file will be altered and checked against Molecule after logic changes '''

import re


class MoleculeTest():
    ''' Object representation of a SMILES code based on symbol by symbol analysis'''

    def __init__(self, smiles, name):

        self.SMILES = self.formatSmiles(smiles)
        self.NAME = name

        self.atomRegex = re.compile(r'[a-zA-Z]')
        self.chargeRegex = re.compile(r'\+|\-')
        self.ATOMS = ['C', 'O', 'N', 'X', 'Z', 'S', 'I',
                      'F', 'c', 'n', 'o', 'x', 'z', 's', 'i', 'f', 'R']
        self.BONDS = ['=', '#']
        self.BRACKETS = ['[', ']']
        self.CHARGES = ['+', '-']
        self.NUMBERS = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
        self.LINEARSYMBOLS = self.ATOMS + self.BONDS + self.BRACKETS + self.CHARGES

        self.RINGS = self.initializeRINGS()
        self.AROMATICINDICES = self.initializeAROMATICINDICES()
        self.CYCLICINDICES = self.initializeCYCLICINDICES()
        self.ringCount = len(self.RINGS)
        self.aromaticCount = 0
        self.nonAromaticCount = 0
        self.determineRingCounts()
        self.RINGDICT = {
            "aromaticRings": self.aromaticCount,
            "nonAromaticRings": self.nonAromaticCount,
            "ringCount": self.ringCount}

        self.SMILES = self.SMILES.upper()
        self.ALCOHOLICINDICES = []
        self.atomData = self.initializeAtomData()
        self.atomCount = len(self.atomData)
        self.bondData = self.initializeBondData()
        self.chargedMol = True if len(
            self.chargeRegex.findall(smiles)) != 0 else False
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

    def initializeRINGS(self):

        openIndex = openPos = -1
        evaluatedNumbers = []
        RINGS = []

        while openPos != len(self.SMILES) - 1:
            openPos += 1
            openSymbol = self.SMILES[openPos]

            if openSymbol in self.ATOMS:
                openIndex += 1
                openAtom = openSymbol

            elif openSymbol == '[':
                openIndex += 1
                chargedGroup = '['

                while self.SMILES[openPos] != ']':
                    openPos += 1
                    chargedGroup += self.SMILES[openPos]
                openAtom = chargedGroup

            if openSymbol in self.NUMBERS and openPos not in evaluatedNumbers:
                closePos = openPos
                closeIndex = openIndex

                while closePos != len(self.SMILES) - 1:
                    closePos += 1
                    closeSymbol = self.SMILES[closePos]

                    if closeSymbol in self.ATOMS:
                        closeIndex += 1
                        closeAtom = closeSymbol

                    elif closeSymbol == '[':
                        closeIndex += 1
                        chargedGroup = '['

                        while self.SMILES[closePos] != ']':
                            closePos += 1
                            chargedGroup += self.SMILES[closePos]
                        closeAtom = chargedGroup

                    if openSymbol == closeSymbol and closePos not in evaluatedNumbers:

                        # RINGPOSITION info in the form [numberPosition, atomIndex, atom]
                        openInfo = [openPos, openIndex, openAtom]
                        closeInfo = [closePos, closeIndex, closeAtom]
                        RINGS.append([openInfo, closeInfo])

                        evaluatedNumbers.append(openPos)
                        evaluatedNumbers.append(closePos)

                        break

        return RINGS

    def numbersHandler(self, pos):

        for ring in self.RINGS:
            openInfo = ring[0]
            closeInfo = ring[1]

            if openInfo[0] == pos:
                return [closeInfo[1], closeInfo[2]]

            elif closeInfo[0] == pos:
                return [openInfo[1], openInfo[2]]

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

    def determineRingCounts(self):

        # Simplification of aromatic/nonaromatic count
        allAtoms = ''.join(self.atomRegex.findall(self.SMILES))

        if allAtoms.islower():
            self.aromaticCount = self.ringCount
            return 0

        elif allAtoms.isupper():
            self.nonAromaticCount = self.ringCount
            return 0

        # Upper is nonaromatic, lower is aromatic
        for ring in self.RINGS:

            # Both are aromatic
            if ring[0][2].islower() and ring[1][2].islower():

                openingNumberPos = ring[0][0]

                if self.SMILES[openingNumberPos+1] in self.ATOMS:

                    if self.SMILES[openingNumberPos+1].islower():
                        self.aromaticCount += 1

                    elif self.SMILES[openingNumberPos+1].isupper():
                        self.nonAromaticCount += 1

                else:
                    self.aromaticCount += 1

            # Both are non aromatic
            elif ring[0][2].isupper() and ring[1][2].isupper():

                openingNumberPos = ring[0][0]

                if self.SMILES[openingNumberPos+1] in self.ATOMS:

                    if self.SMILES[openingNumberPos+1].islower():
                        self.aromaticCount += 1

                    elif self.SMILES[openingNumberPos+1].isupper():
                        self.nonAromaticCount += 1

                else:
                    self.nonAromaticCount += 1

            # They are differnt, must be nonaromatic
            else:
                self.nonAromaticCount += 1

    def initializeAROMATICINDICES(self):

        AROMATICINDICES = []
        atomIndex = -1

        for symbol in self.SMILES:
            if symbol in self.ATOMS:
                atomIndex += 1
                # Lower case atoms are always aromatic
                if symbol.islower() and atomIndex not in AROMATICINDICES:
                    AROMATICINDICES.append(atomIndex)

        return AROMATICINDICES

    def initializeCYCLICINDICES(self):

        evaluatedNumbers = []
        CYCLICINDICES = []

        atomIndex = -1
        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1

            # If a new ring has been found in the SMILES code via a number
            if symbol in self.NUMBERS and pos not in evaluatedNumbers:

                # Index by the scope of paretnehsis groups. I.e. 0th scope is scope where number was found, 1st is an inner parethensis, etc.
                scopeIndices = [[]]
                rBlockCounter = 0
                evaluatedNumbers.append(pos)
                scopeIndices[0].append(atomIndex)
                RNPos = pos + 1
                RNindex = atomIndex
                RN = self.SMILES[RNPos]

                # Symbol is equal to number which opened the ring. Loop untl the same number is encountered by RN
                while RN != symbol:

                    if RN in self.ATOMS:
                        RNindex += 1
                        # If a new nested parenthesis or deeper parenthesis is found, i.e. in between a parenthesis group (...()..,)
                        if len(scopeIndices) == rBlockCounter:
                            # Create a new list tracking the inner indices of that nested parenthesis group
                            scopeIndices.append([RNindex])
                        else:  # If the parenthesis scope is not deeper, then access the current scope and append the atom index to it
                            scopeIndices[rBlockCounter].append(RNindex)

                    if RN == '(':
                        rBlockCounter += 1
                    if RN == ')':
                        del(scopeIndices[rBlockCounter])
                        rBlockCounter -= 1
                    RNPos += 1
                    RN = self.SMILES[RNPos]

                else:
                    evaluatedNumbers.append(RNPos)
                    for scope in scopeIndices:
                        for index in scope:  # Add all indices within that scope to CYCLICindices
                            if index not in CYCLICINDICES:
                                CYCLICINDICES.append(index)

        return CYCLICINDICES

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
        # conversionLegend
        # Br --> X
        # Cl --> Z

        templatePos = -1
        reformattedTemplate = ""

        while templatePos != len(template) - 1:

            templatePos += 1
            if template[templatePos:templatePos+2] == "Br":
                reformattedTemplate += "X"
                templatePos += 1

            elif template[templatePos:templatePos+2] == "Cl":
                reformattedTemplate += "Z"
                templatePos += 1

            else:
                reformattedTemplate += template[templatePos]

        return reformattedTemplate
