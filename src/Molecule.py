""" Represents a SMILES code as its molecular equivalent with atomic properties.

A molecule class contains various containers, holding numerical data, string data, and Atom objects
It derives the bonding schema and cyclic connectivity of the moleucle in space using the SMILES
The class is built on a string decoding algorithm 

Basic usage:

    // For a smiles code
    mol = Molecule('O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O', 'ABEGOH')
    
    // For a functioanl group template
    functionalGroup = Molecule('RC(=O)O', 'Carboxylicacid')

Key Attributes:
    atomData: Dictionary with atomIndex key to Atom Object value pairing
    bondData: Dictionary with atomIndex key to list of Atom Object value pairing
    ringData: Dictionary with ring key to count of specific ring type value

    bondData can be visualized as follows:
    Molecule
    [
        ATOMS→BOND PATHS
        atom1 → [atom2, atom3, atom4…]
        atom2 → [atom1, atom3, atom4…]
        …
        atomN → [atom1, atom2, … atom(N-1)]
    ]

Notes:
    SMILES codes are exported from this module in all upper case
    The module contains Atom objects which are simply symbol index objects 
    SAR = Single Atom Representations
    DLA = Double Lettered Atoms

"""

from Atom import Atom
from helpers import formatSmiles
import re
from collections import OrderedDict


class Molecule():
    """ Object representation of a SMILES code based on symbol by symbol analysis 

    Basic usage:

    Molecule objects:
    >> mol = Molecule('O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O', 'ABEGOH')
    >> print(mol)
    ABEGOH : O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O

    Functional group templats:
    >> functionalGroup = Molecule('RC(=O)O', 'Carboxylicacid')

    This class is utilized for the ifg algorithm
    """

    def __init__(self, smiles, name):
        """ Initialize molecule obejct with bonding and atom data

            Decoding process is as follows:

            1. Decode ring structure of SMILES using numbers and SMILES ring junction pairing
            2. Compute ring counts, aromatic and cyclic properties of atoms
            3. Decode full bonding schema of SMILES using ring pairing and string parsing

            After the decoding process has ran, the Molecule object digitally represents the moleucle
        """

        # Symbol matching lists and regex's for decoding process
        self.ATOM_REGEX = re.compile(r'[a-zA-Z]')
        self.CHARGE_REGEX = re.compile(r'\+|\-')
        self.BOND_REGEX = re.compile(r'\=|\#')
        self.DLA_TO_SAR = {                                  
        "Br": "X",                                       
        "Cl": "Z"
        }
        self.ATOMS =['C', 'O', 'N', 'S', 'I','F', 'P',       # SAR atoms non-aromatic
                    'c', 'n', 'o','s', 'i', 'f', 'p',        # SAR atoms aromatic
                    'R',                                     # R group atom (used for FG templates)
                    'X', 'Z',                                # DLA->SAR converted atoms
                    ]
        self.BONDS = ['=', '#']
        self.BRACKETS = ['[', ']']
        self.CHARGES = ['+', '-']
        self.NUMBERS = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
        self.LINEARSYMBOLS = self.ATOMS + self.BONDS + self.BRACKETS + self.CHARGES

        # Input data
        self.SMILES = formatSmiles(smiles)
        self.NAME = name

        # Ring Containers
        self.RING_OPEN_POSITIONS = []
        self.RING_SELF = {}
        self.RING_CLOSE_POSITIONS = []
        self.RING_COMPLEMENTS = {}
        self._RING()                # Initalize Ring Containers

        # Ring Index Containers
        self.AROMATICINDICES = []
        self.CYCLICINDICES = []
        self._INCIDES()             # Initalize Ring Index Containers

        # Collect ring data from ring containers
        self.ringCount = len(self.RING_SELF) / 2
        self.aromaticCount = 0
        self.nonAromaticCount = 0
        self._RING_COUNTS()         # Compute Ring counts

        # Dictionary mapping count by name to amount
        self.ringData = {
            "aromaticRingCount": self.aromaticCount,
            "nonAromaticRingCount": self.nonAromaticCount,
            "ringCount": self.ringCount
        }

        # Atom and Bonding data
        self.ALCOHOLICINDICES = []                  # Holds indices of oxygens with alcoholic properties
        self.SMILES = self.SMILES.upper()           # After aromaticity analysis, convert SMILES to full uppercase for simplicity
        self.atomData = self.initializeAtomData()   # Atom index to symbol pairing dictionary
        self.atomCount = len(self.atomData)         
        self.bondData = self.initializeBondData()   # Atom index to bonded atoms pairing dictionary
        self.chargedMol = (
                True if len(self.CHARGE_REGEX.findall(smiles)) != 0 
                else False
            )
        self.AMINOACID = (
                True if len(re.compile(r'\[[nN]H[23]?\+\]').findall(smiles)) != 0 
                else False
            )

        self.HALOGENS = {
            "Bromine" : len([atom for atom in self.atomData.values() if atom.symbol == "X"]),
            "Chlorine" : len([atom for atom in self.atomData.values() if atom.symbol == "Z"]),
            "Iodine" : len([atom for atom in self.atomData.values() if atom.symbol == "I"]),
            "Fluorine" : len([atom for atom in self.atomData.values() if atom.symbol == "F"])
        }

        self.ATOMCOUNTS = {
            "Br" : len([atom for atom in self.atomData.values() if atom.symbol == "X"]),
            "Cl" : len([atom for atom in self.atomData.values() if atom.symbol == "Z"]),
            "I" : len([atom for atom in self.atomData.values() if atom.symbol == "I"]),
            "F" : len([atom for atom in self.atomData.values() if atom.symbol == "F"]),
            "C" : len([atom for atom in self.atomData.values() if atom.symbol == "C"]),
            "O" : len([atom for atom in self.atomData.values() if atom.symbol == "O"]),
            "N" : len([atom for atom in self.atomData.values() if atom.symbol == "N"]),
            "S" : len([atom for atom in self.atomData.values() if atom.symbol == "S"]),
            "P" : len([atom for atom in self.atomData.values() if atom.symbol == "P"]),
        }

    def __str__(self):
        """ String representation of a Molecule SMILES object 
        """
        
        return self.NAME + " : " + self.SMILES

    def getSymbolDict(self):
        """ Create a dictionary of atom index to atom symbol of the molecule's atoms based on its current state
        """

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
        atomIndex = -1                                          # 0 based indexing of atoms

        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:                            # Exclusively analyze atoms
                atomIndex += 1
                atomSymbol = symbol
                
                if self.SMILES[pos-1] == '[':                   # If atom charged, collect its bracketed group
                    atomSymbol = self.getChargedGroup(pos-1)    # Use opening bracket for charge group

                atom = Atom(atomIndex, atomSymbol)              # Create new atom objects
                atomData[atomIndex] = atom                      # Atom index to atom object pairing
                
                if symbol == 'O':                               # Alochols stem from oxygens
                    self.determineAlcoholGroup(pos, atomIndex)

        return atomData

    def initializeBondData(self):
        """ Returns the bonding schema for a given SMILES

            Psuedo-code
            For each atom in the SMILES
                Determine the left bonds
                Determine the right bonds
                Combine both bond pathways into a single pathway that stems from current atom
                Assign the atomIndex of bondData to its bonding pathways
        """

        bondData = {}
        atomIndex = -1                                  # 0 indexed atom schema

        for pos, symbol in enumerate(self.SMILES):

            if symbol in self.ATOMS:
                atomIndex += 1
                bonds = []                              # Bonds which stem from a single atom

                if self.SMILES[pos-1] == '[':           # Charged atom case

                    leftBond = self.getLeftBond(        
                            LNPos=pos-1,                # Start left pathway at open bracket
                            LNIndex=atomIndex
                    )

                    rightBonds = self.getRightBonds(   
                            RNPos=pos+2,                # Start right pathway at closed bracket
                            RNIndex=atomIndex
                    )

                else:                                   # Uncharged atom case
                    leftBond = self.getLeftBond(        
                            LNPos=pos, 
                            LNIndex=atomIndex
                    )

                    rightBonds = self.getRightBonds(
                            RNPos=pos, 
                            RNIndex=atomIndex
                    )

                if leftBond:                            # Left bond is always one bond if it exists
                    bonds.append(leftBond)

                if rightBonds:                          # Right bonds can be more than 1 (parenthesis and number cases)
                    for rightBond in rightBonds:
                        bonds.append(rightBond)

                bondData[atomIndex] = bonds             # Assign the atom index to its bonded neighbors with the bondData dictionary
                del(bonds)

        return bondData

    def getLeftBond(self, LNPos, LNIndex):
        """ Return the single Atom object that is left bonded with respect to a specific string position
            Return an empty array [] if no LN exists from a specific string position

            LNPos (int) : string position of the atom whose left bond is to be found
            LNIndex (int) : index position of that atom within the smiles code

            Notes:
                LN stands for Left Neighbor
        """

        LNPos -= 1                  # Begin one position to left of given atom position
        if LNPos < 0:               # Outside scope of SMILES means no left bonds
            return []
        LN = self.SMILES[LNPos]     # Initalize LN as the first symbol to left of positon given

        explicitBond = (            # Explicit double/triple left bond case
            LN 
            if LN in self.BONDS     # Save results of bond if LN is initially a bond
            else ""
        )
        if explicitBond:            # If an explicit bond is found
            LNPos -= 1              # Move to next symbol, bond is saved in explicitBond
            LN = self.SMILES[LNPos]

        scope = 0                   # Parenthesis depth counter is 0 index based
        leftBond = []               # Left bond list

        
        while LN not in self.LINEARSYMBOLS or scope > 0:    # Loop until a LINEARSYMBOL is found on the 0th scope (i.e. same scope as inital atom)

            if LN == ')' and scope == -1:                   # Double parenthesis case i.e X(Y)(Z) starting within second parenthesised Z
                scope = 0                                   # Reset scope to allow second atom Z to locate X as its proper left-bonded neighbor
            if LN == ')':                                   # Increment over deeper parenthesis groups
                scope += 1
            elif LN == '(':                                 # Decrement into higher parenthesis gropus
                scope -= 1
            elif LN in self.ATOMS:                          # Decrement atomic index for each atom found
                LNIndex -= 1
            LNPos -= 1                                      # Decrement to next left symbol
            if LNPos < 0:                                   # Outside scope of SMILES means no left bonds
                return []
            LN = self.SMILES[LNPos]                         # Update LN

        else:                                               # LN in this scope is the correct left-bonded atom of given arguments
            LNIndex -= 1                                    # Decrement to find corrected LN index
            if LN == ']':                                   # Charged LN case
                LN = self.getChargedGroup(LNPos)            # Find LN directly from final position 

            leftBond = Atom(LNIndex, explicitBond + LN)     # Create new bonded atom at LNindex with the explicitBond and the determined LN

        return leftBond                                     

    def getRightBonds(self, RNPos, RNIndex):
        """ Returns an atom list of all right bonds which stems from a specific string position
            Returns an empty list if no right bond exists from the given position


            RNPos (int) : The position of the atom in the SMILES code whose right bonds are to be found
            RNIndex (int) : The index position of the atom in the SMILES code

            Notes:
                This is a recursive function
                RN stands for Right Neighbor
                Right groups split to different potions of the smiles code (dynamic symbols)
                The two dynamic symbols are () and numbers. Both are handled accordingly here. 
        """

        RNPos += 1                       # Begin at one position to right
        if RNPos >= len(self.SMILES):    # Outside SMILES scope means no right bond
            return []
        RN = self.SMILES[RNPos]          # Initalize RN as the symbol one position to right of positional arguments

        if RN == ')':                    # A closing parenthesis explicitly means no right bonds
            return []

        scope = 0                        # Parenthesis depth counter is 0 index based
        rightBonds = []                  # Right bonding list

        while RN not in self.LINEARSYMBOLS or scope > 0:    # Loop until a LINEARSYMBOL is found on the 0th scope (i.e. same scope as inital atom)

            if RN == '(':                                   # Nested parenthesis case

                innerParenthBond = self.getRightBonds(      # Only allow single depth right bond neighbors (i.e. first bond on 1st scope)
                    RNPos, RNIndex
                ) if scope == 0 else ''                     # Nonzero scope implies deeper RN's unconnected to the initial posiiton relative to arguments

                if innerParenthBond:                        # Variable is a list of bonds, but only first bond is valid instead of deeper scoped bonds
                    rightBonds.append(innerParenthBond[0])

                scope += 1                                  # Increment scope depth

            if RN == ')':                                   # Decrement scope for right bond upon closing parenthesis
                scope -= 1
            if RN in self.ATOMS:                            # Increment atom index each atom encountered
                RNIndex += 1
            if scope == 0 and RN in self.NUMBERS:           # Every number always has a bonded partner
                numGroup = self.numbersHandler(RNPos)       # Handle ring junction bond via numbersHandler
                rightBonds.append(numGroup)                 # Add bonded group to list of bonds

            RNPos += 1                                      # Increment position to process next symbol
            if RNPos >= len(self.SMILES) or scope < 0:      # Outside of scope 
                return rightBonds                           # Return current state of bonded list
            RN = self.SMILES[RNPos]                         # Update RN

        else:                                               # RN is the final right bonded neighbor within this scope
            RNIndex += 1                                    # Increment atomic index to correct index offset

            if RN in self.ATOMS:                            # Generic atom right bond case
                rightBonds.append(Atom(RNIndex, RN))

            elif RN in self.BONDS:                          # Explicit right bond case
                explicitBond = RN                           # RN must be explicit bond
                RNPos += 1
                RN = self.SMILES[RNPos] 

                if self.SMILES[RNPos] == '[':               # Explicit bond to charge group case
                    RN = self.getChargedGroup(RNPos)

                rightBonds.append(Atom(RNIndex, explicitBond + RN))     # Create new atom bond with explicit bond and RN

            elif RN == '[':                                             # Explicit charge group right bond case
                chargedGroup = self.getChargedGroup(RNPos)
                rightBonds.append(Atom(RNIndex, chargedGroup))          # Create new atom with full charged group no explicit bond case

        return rightBonds

    def getChargedGroup(self, pos):
        """ Returns a charged group based on a selected bracket position in the its smiles code
            pos (int): string position of an openeing or closing bracket in the smiles code
        """

        chargedGroup = ''                           # Resultant charged group string

        if self.SMILES[pos] == ']':                 # Closing (Reversed) case
            while self.SMILES[pos] != '[':
                chargedGroup += self.SMILES[pos]
                pos -= 1
            return '[' + chargedGroup[::-1]         # Reverse resultant string

        if self.SMILES[pos] == '[':                 # Opening (Forwards) case
            while self.SMILES[pos] != ']':
                chargedGroup += self.SMILES[pos]
                pos += 1
            return chargedGroup + ']'

        if self.SMILES[pos] not in self.BRACKETS:   # Capture errors in positions
            raise ValueError(
                "The position passed does not point to a bracket in the smiles code"
            )

    def _RING(self):
        """ Initalizes the four ring data containers
            RING_SELF (dictionary): string position of number to its direct atom pair (i.e. open number position gives opening atom data)
            RING_COMPLEMENTS (dictionary): string position of a number to its complementary atom pair (i.e. open number position gives closing atom data)
            RING_OPEN_POSITIONS (list): string positions of the numbers which opened a junction
            RING_CLOSE_POSITIONS (list): string positions of the numbers which closed a junction


            Notes: 
                Ring junctions open with an arbitrary number and close with the same number
                The atoms to the left of each number are at the opening/closing atoms at ring junction
        """

        atomIndex = -1              # 0 based indexing of atoms
        evaluatedNumbers = {}       # Number to string position in smiles code
        atomSymbol = ''

        for pos, symbol in enumerate(self.SMILES):              # Loop over SMILES 

            if symbol in self.ATOMS:
                atomIndex += 1
                atomSymbol = symbol

            if symbol in self.NUMBERS:                          # Symbol is a number in this scope

                if self.SMILES[pos-1] == ']':                   # Charged group at ring junction
                    atomSymbol = self.getChargedGroup(pos-1)    # Retrieve charged group from closing bracket

                atom = Atom(atomIndex, atomSymbol)              # Ring junction atom object
                self.RING_SELF[pos] = atom                      # Opening ring junction string position to opening atom pair

                if symbol in evaluatedNumbers:                  # If a closing ring junction has been located

                    initalNumberPos = evaluatedNumbers[symbol]                  # Inital string position of number where ring junction opened
                    self.RING_COMPLEMENTS[initalNumberPos] = atom               # Opening ring junction string position to closing atom pair
                    self.RING_COMPLEMENTS[pos] = self.RING_SELF[initalNumberPos]# Closing ring junction string position to opening atom pair
                    self.RING_CLOSE_POSITIONS.append(pos)                       # Position of where the ring ends in the SMILES
                    del(evaluatedNumbers[symbol])                               # Close ring path once completed to allow other paths with the same number
                    continue                                                    # Process the next ring

                self.RING_OPEN_POSITIONS.append(pos)            # Position of where the ring starts in the SMILES
                evaluatedNumbers[symbol] = pos                  # Number to positioning key:value pairing

    def numbersHandler(self, pos):
        """ Return the partner of a given position to retrieve what atom a specific number is connected to

            pos (int) : Position of a ring whose partner atom is to be determined
        """

        return self.RING_COMPLEMENTS[pos]

    def determineAlcoholGroup(self, pos, atomIndex):
        """ Validates if an given string position and atom index in a Molecule correspond to an alochol group
            Valid alcoholic oxygen indicies will be added to the ALOCHOLICINDICIES list 

            pos (int) : position of the oxygen to be checked for alocholic properties
            atomIndex (int) : index position of the oxygen

            Notes:
                Only four cases where an alcohol is determined. 
                Currently alochol cases are hard coded here
                The index of the oxygen is recorded as being alocholic using the ALOCHOLICINDICES list
        """

        if pos == len(self.SMILES) - 1:                         # Final atom alcohol group case
            if(
                self.SMILES[pos-1].upper() == 'C'               # Can be connected to an linear carbon
                or self.SMILES[pos-1] in self.NUMBERS           # Or connected at a ring junction
                or self.SMILES[pos-1] == ')'                    # Or connected at a normal junction
            ):
                self.ALCOHOLICINDICES.append(atomIndex)         

        elif(                                                   # First atom alcohol group case
            atomIndex == 0                                      # At index 0
            and self.SMILES[pos+1].upper() == 'C'               # Can only be connected to a linear carbon
        ):
            self.ALCOHOLICINDICES.append(atomIndex)            

        elif self.SMILES[pos-1:pos+2] == '(O)':                 # Lone alcohol group case
            self.ALCOHOLICINDICES.append(atomIndex)

        elif(                                                   # Ending parenthesis alcohol group
            self.SMILES[pos+1] == ')'                           # Must be directly next to closing parenthesis
            and self.SMILES[pos-1] not in self.BONDS            # Not connected to a double/triple bond
            and self.SMILES[pos-1] not in self.BRACKETS         # And not connected to a charge group
        ):
            self.ALCOHOLICINDICES.append(atomIndex)

        return 

    def _RING_COUNTS(self):
        """ Loops over Ring containers to determine the counts of aromatic and non aromatic rings in the molecule

            Notes:
                Aromaticity takes priority (i.e. if a single atom within ring is aromatic, take ring as being aromatic)
                Uses the two atoms involved in the number ring positioning to detemine what kind of ring the atoms are apart of
                Some extra neighbor checking is required to confirm aromatic/non aromatic presence fully
        """

        allAtoms = ''.join(self.ATOM_REGEX.findall(self.SMILES))    # Simplification of aromatic/nonaromatic count

        if allAtoms.islower():                                      # All lower case atoms means all counted rings are aromatic
            self.aromaticCount = self.ringCount
            return 

        elif allAtoms.isupper():                                    # All upper case atoms means all counted rings are non-aromatic
            self.nonAromaticCount = self.ringCount
            return 

        for pos in self.RING_OPEN_POSITIONS:                        # Loop over the opening ring junction string positions in a Molecule
            ringOpen = self.RING_SELF[pos]                          # Get the opening atom data
            ringClose = self.RING_COMPLEMENTS[pos]                  # Get the closing atom data

            if(                                                     # Both atomic junctions are aromatic case
                ringOpen.symbol.islower() 
                and ringClose.symbol.islower()
            ):

                if self.SMILES[pos+1] in self.ATOMS:                # Check if first symbol deeper in ring is an atom

                    if self.SMILES[pos+1].islower():                # Another lower atom in ring confirms aromatic ring
                        self.aromaticCount += 1

                    elif self.SMILES[pos+1].isupper():              # Non-aromatic atom in ring means polycyclic non-aromatic and aromatic ring interaction
                        self.nonAromaticCount += 1

                else:                                               # If next symbol is not an atom, then the ring is aromatic
                    self.aromaticCount += 1                         

            elif(                                                   # Both atomic junctions are non aromatic case
                ringOpen.symbol.isupper() 
                and ringClose.symbol.isupper()
            ):

                if self.SMILES[pos+1] in self.ATOMS:                # Check if first symbol deeper in ring is an atom

                    if self.SMILES[pos+1].islower():                # Another lower atom in ring confirms aromaic ring
                        self.aromaticCount += 1

                    elif self.SMILES[pos+1].isupper():              # Otherwise another upper atom confirms non-aromatic ring
                        self.nonAromaticCount += 1

                else:                                               # Ring is non-aromatic otherwise
                    self.nonAromaticCount += 1                      

            else:                                                   # Different ring junction atom types means take non-aromatic priority
                self.nonAromaticCount += 1

    def _INCIDES(self):
        """ Determines the cyclic/aromatic atom indicies inside the SMILES

            AROMATICINDICES (list): List of atom indicies which are apart of an aromatic ring
            CYCLICINDICES (list): List of atom indicies which are apart of a ring

            Notes:
                Scope is the parenthesis group level, counted relative to numbered opening position, which starts at 0th index
                If a ring ceases inside a parenthesis group, then all indicies in that scope and above it are contained in the ring
                Therefore, tracking of indices on specific scope levels (parenthesis scopes) is necessary to obtain the most accurate
                ring indices

                If a number ends within a depth deeper than 0, i.e. ...1... (... (...1)...),
                Then all all indices in in between must be cyclic
                Some nested scopes may go deeper, but may not conclude the ring.
                For example, ...1... (... (...) ... 1). The middle parenthesis is not apart of the ring structure related to the number 1.
                However, the number 1 still closes within a nested parenthesis. Because the depth of conclusion is not known before this algorithm is run
                All depths within indivudal depth counts must be tracked to obtain accurate cyclic index information

                Form of scopeIndices
                scopeIndicies =
                [
                    [index1,index2...],        scope = 0           (same scope as where the number appears)
                    [index3,index4...],        scope = 1           (one scope deeper from where number appears)
                    ...
                    [indexN...]                scope = N-1         (N'th scope deeper from where number appears)
                ]
        """

        atomIndex = -1                  # 0 based indexing of SMILES
        evaluatedNumPositions = []      # String positions of numbers in SMILES (opening/closing ring junctions)

        for pos, symbol in enumerate(self.SMILES):              # Loop over SMILES

            if symbol in self.ATOMS:
                atomIndex += 1

            if symbol.islower() and atomIndex not in self.AROMATICINDICES:       # Lower case letters are aromatic atoms
                self.AROMATICINDICES.append(atomIndex)
                continue

            
            if symbol in self.NUMBERS and pos not in evaluatedNumPositions:      # If a new ring has been found in the SMILES code via a number

                evaluatedNumPositions.append(pos)               # Current position is a number
                scopeIndices = [[]]                             # ScopeIndicies list for cyclic index tracking
                scope = 0                                       # 0 index scope for number scope
                scopeIndices[0].append(atomIndex)               # Add first index to 0th scope level
                RNPos = pos + 1                                 # Start RN next to opening ring number
                RNindex = atomIndex                             # Start RN index at most recent atom index
                RN = self.SMILES[RNPos]                         # Get symbol from SMILES
               
                while RN != symbol:                             # Loop until ring closes (symbol is opening number within this scope)

                    if RN in self.ATOMS:                        
                        RNindex += 1                            # Increment atom index

                        if len(scopeIndices) == scope:          # If this atom is apart of a deeper parenthesis scope, i.e. 1(...
                            scopeIndices.append([RNindex])      # Then create a new list for tracking the indices part that scoped group
                        else:                                   # Otherwise, the atom is still apart of the same scope
                            scopeIndices[scope].append(RNindex) # Add atom index to scoped list

                    if RN == '(':                               # Open parenthesis increments scope depth 
                        scope += 1                  

                    if RN == ')':                               # Close parenthesis decrements scope depth and removes deepest scope from ring 
                        del(scopeIndices[scope])               
                        scope -= 1                              

                    RNPos += 1                                  # Go to next symbol position
                    RN = self.SMILES[RNPos]                     # Get next symbol

                else:                                           # Closing ring number has been located
                    evaluatedNumPositions.append(RNPos)         # Add number position to evaluated number positions
                    for scope in scopeIndices:                  # Every unique index in scopeIndicies is cyclic, AKA apart of a ring
                        for index in scope:                     
                            if index not in self.CYCLICINDICES: 
                                self.CYCLICINDICES.append(index)

