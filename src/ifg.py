""" Algorithm which analyzes the atom and bond data of a given SMILES code to produce counts of pre-defined organic functional groups.

    This algorithm is built off of the Molecule object, which has two important containers: atomData and bondData. These two
    dictionaries possess the data of how the molecule atoms are connected. bondData contains all the possible bond paths that stem from
    every atom in the SMILES code. The same is true of decoded functional group templates from the Molecule class.
    Therefore, if a functional group template can model itself upon at least one of those paths in the SMILES code, 
    then that functional group must be present inside of the SMILES code. The concept of a valid path for a functional group
    within the SMILES code string is the solution to this algorithm.

Key Attributes:

    The following two functional group containers have keys for the counts of their that specific SMILES code functional groups, 
    including extra classification like cyclic or aromatic. Also includes ring data and alcohol counts in both. They differ slightly:


    allFgs (Dict) : The dictionary which counts all functional groups.
    "All" is defined as allowing overlapping functional groups, such as a ketone inside of an ester. 
    
    preciseFgs (Dict) : The dictionary which counts precise functional groups.
    "Precise" is defined as allowing overlapping functional groups. This means ketones would not be seen inside esters, 
    nor tertiary amines inside amides, and so on. Only the largest functional groups will appear within this dictionary
"""


from Molecule import Molecule
import re
import os


class ifg(Molecule):
    """ The class algorithm which identifies the functional group counts of a SMILES code. """

    def __init__(self, SMILES, REFCODE):
        """ Initalizes SMILES code Molecule representation, then processes the functional group algorithm on that molecule. 

            SMILES (string) : A valid simplified molecular input line entry system (SMILES) code
            REFCODE (string) : The referenced code for this particular SMILES code
        
            Process of algorithm:
            
            1. Decode the input SMILES code into a digital Molecule
            2. Determine the functional groups from the digital Molecule

            Notes:
                invoking "self." within this class after the super().__init__(SMILES, REFCODE) function call will allow 
                access to all the molecular properties for the specific input SMILES codes.
        """

        super().__init__(SMILES, REFCODE)                               # Create Molecule object from input SMILES code
        self.allFgs = self.findFunctionalGroups()                       # Entry point of full IFG algorithm
        self.preciseFgs = self.findPreciseGroups(self.allFgs)           # Filter output to obtain precise group counts

    def findFunctionalGroups(self):
        """ Return a list of Molecule objects which represent the functional groups of a particular SMILES code

            Notes:
                FGs as Molecules inside of the SMILES allows the atomic indicies to
                reflect the positioning within the SMILES
            """

        functionalGroups = []

        for atom in self.atomData.values():                     # Loop over all atoms in the SMILES
            groups = self.whichGroup(atom)                      # Determine the functional groups which stem from an atom

            if atom.symbol == 'N':                              # Special handling for nitrogens
                primaryAmine = self.detetminePrimaryAmine(atom.index)

                if primaryAmine:
                    functionalGroups.append(primaryAmine)

            for group in groups:                                # For each group (each are Molecule objects) found to stem from an atom
                functionalGroups.append(group)                  # Add it to the full list of identified FG

        self.repetitionFilter(functionalGroups)                 # Remove repeated groups
        self.priorityFilter(functionalGroups)                   # Remove overlapped groups based on priority of R groups 

        self.determineCyclicGroups(functionalGroups)            # Identify the groups which are cyclic by index analysis
        self.determineAlcoholGroups(functionalGroups)           # Create Molecule objects for identified alcohols

        return functionalGroups

    def whichGroup(self, atom):
        """ Returns the list of functional groups as Molecule objects that stem from a given atom in the SMILES code

            atom (type Atom) : An atom object who is to be analyzed for being contained in a possible functional group

            Process is as follows:
                Filters functional groups which do not contain given atom symbol
                Selects an expansion point from within a potential functional group to start expansion from 
                    (i.e. Which C to choose in [R]C(=O)NC(=O)[R], 1st or 2nd)
                Creates Molecule object of functional group
                Calls expandgroup to complete path validation 
            
            Notes:
                Loops over every possible functional group in FGlist
                Every atom has a possibility to be branched into a functional group. Therefore, every atom must be looked at individually.
                R group paths are elimated because they cannot be expanded from to validate the functional group. Only main group atoms can be expanded from
                Function currently ONLY SUPPORTS SINGLE POINT ALCOHOL EXPANSION (i.e. only 1 alcohol is allowed in FG)
        """

        matches = []
        for line in open(os.getcwd() + '/src/resources/FGlist.txt', 'r'):           # Loop over every functional group
            lineInfo = re.compile(r'\S+').findall(line)
            lineInfo[0] = lineInfo[0].replace('[R]', 'R')
            FGtemplate = Molecule(lineInfo[0], lineInfo[1])                         # FG Molecule object

            if(                                                                     # Non-charged SMILES cannot contain a charged FG
                len(self.CHARGE_REGEX.findall(atom.symbol)) == 0                    # Charged FG's
                and len(self.CHARGE_REGEX.findall(FGtemplate.SMILES)) != 0          # In a non-charged SMILES
            ):
                continue                                                            # Are skipped

            if(                                                                     # Alcoholic FG's can only stem from alcoholic oxygens
                FGtemplate.ALCOHOLICINDICES                                         # Alcoholic FG
                and atom.index not in self.ALCOHOLICINDICES                         # With a Non-alcoholic oxygen atom
            ):
                continue                                                            # Are skipped


            if atom.symbol in FGtemplate.SMILES:                                    # If the atom symbol argument is inside the FGtemplate
                expansionPoint = 0                                                  # String position within FG where the atom is located

                for tempAtom in FGtemplate.atomData.values():                       # R group pathing is elimated in FGtemplate
                    if tempAtom.symbol == 'R':                                      # Locate R groups in FGtemplate
                        FGtemplate.bondData[tempAtom.index].clear()                 # Remove the bonding paths from the R groups

                for tempAtom in FGtemplate.atomData.values():                       # Find expansion indicies or points inside FGtemplate

                    if not FGtemplate.ALCOHOLICINDICES:                             # Non-alcoholic FG's can expand from any main group atom
                        if atom.symbol == tempAtom.symbol:                          # When the given atom is located within the FG
                            expansionPoint = tempAtom.index                         # Set the string expansion point at that atom within the FG
                            break
                    else:                                                           # Alcoholic FG's can only expand from its alcohol group
                        expansionPoint = FGtemplate.ALCOHOLICINDICES[0]             # Only one alcohol is considered inside of the FG
                        break

                                                                                    # Set FGtemplate for expandGroup
                FGtemplate.atomData[expansionPoint].index = atom.index              # Match the FG atom index with its SMILES code index equivalent
                self.removeBond(expansionPoint, FGtemplate)                         # Remove the starting atom from all bonding paths, as it will be stemmed from

                expandGroup = self.expandGroup( atom,                               # Attempt to expand the FGtemplate inside of the SMILES code
                                                expansionPoint, 
                                                FGtemplate, 
                                                None, 
                                                atom.index
                )
                if expandGroup:                                                     # If expansion was successful
                    matches.append(FGtemplate)                                      # Add the Molecule group, now with corresponding indicies to the SMILES code, to the lits of matches

        return matches

    def expandGroup(self, atom, expansionPoint, template, skipIndex=None, smilesIndexInit=None):
        """ Return true or false for if a specifed atom can be expanded into the passed in template.

            atom (type Atom) : the SMILES code atom which is to be analyzed for having a valid bond with respect to the template
            expansionPoint (int) : The position inside of the template for which atom is to be decoyed as, or tested aginst for bonds matching inside of the template. 
            template (type Molecule) : The functional group for which the given atom will be checked against.
            skipIndex (int, default=None) : A template bond for which to skip over
            smilesIndexInit (int, default=None) : The indices confirmed in the SMILES code to be apart of a functional group, with designated positions

            Notes:
                This is a recusrive algorithm. It takes two atoms, which have identical symbol, and determines if a bonded path is equivalent in both. 
                This algorithm answers the question: Is this atom, in the smiles code, bonded in the same that another atom, in the temaplte, is bonded?
                If an atom in the smiles code can trace its own bonds to mirror those inside of a functional group, then that functional must exist within the smiles code

                The recursion comes from multiple atoms needed to be called upon. Example of how the algorithm will work

                    Assume a nitrogen has been found, and we would like to validate if that nitrogen is part of an amide structure
                    template: [R]C(=O)N([R])[R]
                    atom: Atom(symbol=N, index=5)

                    In whichGroup, this N is prepared to "fit" the template. expandGroup will retrieve an nitrogen atom object. The atom index, 5 in this case,
                    will be used to retrieve its bonded content from the its smiles molecule, which in this case is accessed by self (remember super call in __init__ method)

                    The template has its own bond schema as well. Given an expansion point, i.e. where the atom of interest in the TEMPLATE is located, 
                    the necessary bonds to complete a functional group can be retirved with the temaplte.bondData attribute. 

                    With both pieces of bondData from each atom, it is now possible to check if valid bonding is present inside of the smiles code 
                    for a given functional group, represented in template

                    All template bonds are looped over. The skipIndex inside of template is skipped.
                    This is implemented to prevent a template bond from being validated by the bond it was just stemmed from during a recursive call

                    The specific bond in the template is now checked against the available smiles bonds. If a unique and unsused smiles bond (i.e smiles atom)
                    is available for usage by the template (i.e. An availble carbon was just found for the nitrogen). 
                    Then the bond path stemming from that atom must also be validated

                    Continually, this algorithm first validates all of the main group atoms which are necessary for a functional group are indeed paired
                    inside of the smiles code. Only after are the RGroups validated. In this way, both types of atoms are seperated, and main/Rgroup atoms get
                    distinct delegation. There is no extra embedded crossover, they are handled independently for each atom in the template, with respect to the SMILES

                    Remember, if a path is validated and a new path has begun, the bond representing the atom from which we just came from should not be used
                    It already has a place in the template, and is being represented in the smiles code. Therefore, one path is validted by virtue of being called upon

                    If all R groups and previous paths are valid, then expand true finally returns true: This means that the entire bond path from a specific atom in a tempalte
                    can be mirrored into the smiles code

                    If at any point expand group fails and returns false, it means that a specifc atom could not be branched into the functional group via its smiles bonding schema.

                    For all valid atoms found, the template atomData is overwritten with proper SMILES index representations. Therefore, by appending template in whichGroup, 
                    a functional group with indices pointing to the atoms within its SMILES code is saved and collected. 
        """

        smilesIndices = []  # Smiles indicies which have already been used in template
        smilesIndices.append(smilesIndexInit)
        templateBonds = template.bondData[expansionPoint]
        mainTemplateBonds = self.filterRgroup(templateBonds)
        smilesBonds = self.bondData[atom.index]
        templateAtoms = template.atomData

        # Not enough bonds at desposal of smiles code to satisfy template means False match
        if len(templateBonds) > len(smilesBonds):
            return False

        # Main group Template Analysis
        for tempBond in mainTemplateBonds:
            tempIndex = tempBond.index
            if tempIndex == skipIndex:
                continue
            tempSymbol = tempBond.symbol
            for smilesBond in smilesBonds:

                smilesIndex = smilesBond.index
                smilesSymbol = smilesBond.symbol

                if tempSymbol == smilesSymbol and smilesIndex not in smilesIndices:

                    path = self.expandGroup(
                        smilesBond, tempIndex, template, expansionPoint, atom.index)

                    if path:
                        smilesIndices.append(smilesIndex)
                        templateAtoms[tempIndex].index = smilesIndex
                        break
                    else:
                        continue
            else:
                return False

        # Rgroup Template Analysis
        else:

            RgroupCounter = -1
            Rgroups = self.getRgroups(templateBonds)

            if not Rgroups:
                return True

            if len(smilesBonds) > len(mainTemplateBonds):

                for smilesBond in smilesBonds:

                    smilesIndex = smilesBond.index
                    smilesSymbol = smilesBond.symbol

                    if smilesIndex in smilesIndices or smilesSymbol[0] in self.BONDS:
                        continue
                    RgroupCounter += 1
                    if RgroupCounter < len(Rgroups) - 1:
                        templateAtoms[
                            Rgroups[RgroupCounter].index
                        ].index = smilesIndex
                    else:
                        templateAtoms[
                            Rgroups[RgroupCounter].index
                        ].index = smilesIndex
                        break
                else:
                    return False

                return True

            else:
                return False

    def getRgroups(self, atomSet):
        """ Returns a dictionary or list of all R group atoms in a given atomSet.

            atomSet (list || dict) : Data container which has Atom objects, either as values or in a list
        """

        if isinstance(atomSet, list):
            return [atom for atom in atomSet if atom.symbol == 'R']
        if isinstance(atomSet, dict):
            return {atom.index: atom for atom in atomSet.values() if atom.symbol == 'R'}

    def filterRgroup(self, atomSet):
        """ Returns a dictionary or list of all non R group atoms in a given atomSet.

            atomSet (list || dict) : Data container which has Atom objects, either as values or in a list
        """

        if isinstance(atomSet, list):
            return [atom for atom in atomSet if atom.symbol != 'R']
        if isinstance(atomSet, dict):
            return {atom.index: atom for atom in atomSet.values() if atom.symbol != 'R'}

    def removeBond(self, index, template):
        """ Removes a bond path from a template.

            index (int): index of atom in template whose bonds are to be removed
            template (type Molecule) : The functional group template who contains the bond to be removed

            Notes:
                When other atoms in the template look at their bonds after this function call, 
                they will no longer see the atom whose index is the index argument passed in
        """

        for i, bonds in enumerate(template.bondData.values()):
            for j, bond in enumerate(bonds):
                if bond.index == index:
                    del(template.bondData[i][j])
        return 

    def repetitionFilter(self, functionalGroups):
        """ Returns a filtered list of functionalGroups for repeated groups

            functionalGroups (list) : List of Molecule obejects that represent the functional groups in a smiles code

            Notes:
                An FG is identified as repeated if all the atoms within each one are identical
                If every atom is identical in both functional groups and they are unique FGs identified individually, then they are repeated
        """

        index = -1                                      

        while index != len(functionalGroups) - 1:       

            index += 1
            group = functionalGroups[index]
            groupAtoms = group.atomData.values()

            for compareIndex, compareGroup in enumerate(functionalGroups):

                compareAtoms = compareGroup.atomData.values()

                if(
                    all(i in compareAtoms for i in groupAtoms) 
                    and all(i in groupAtoms for i in compareAtoms) 
                    and compareIndex != index
                ):
                    del(functionalGroups[compareIndex])
                    index = - 1
                    break

        return 

    def determineCyclicGroups(self, functionalGroups):
        """ Label functional groups with ring classification based on the ring structure of molecule (self) 

            functionalGroups (list) : List of Molecule obejects that represent the functional groups in a smiles code

            Notes: 
                Some rings are directly connected to their opposite type, i.e aromatic-nonAromatic
                To be the most accurate, a tally system was created to tally how many valid bonds fall under a certain classification.
        """

        for group in functionalGroups:
            groupAtoms = group.atomData
            aromaticCount = cyclicCount = 0
            for templateIndex, smilesAtom in groupAtoms.items():

                if smilesAtom.symbol not in self.LINEARSYMBOLS:
                    continue

                # If maingroup atom does not have any cyclic property, then skip it
                if smilesAtom.index not in self.AROMATICINDICES and smilesAtom.index not in self.CYCLICINDICES:
                    continue

                # If two consecutive atoms, bonded together, are part of a ring, add it to the "count"
                for templateBond in group.bondData[templateIndex]:
                    smilesBond = groupAtoms[templateBond.index]
                    indicies = [smilesAtom.index, smilesBond.index]
                    if all(i in self.AROMATICINDICES for i in indicies):
                        aromaticCount += 1
                    elif all(i in self.CYCLICINDICES for i in indicies):
                        cyclicCount += 1
            else:
                if aromaticCount or cyclicCount:
                    if aromaticCount >= cyclicCount:
                        group.NAME = 'Aromatic' + group.NAME
                    elif aromaticCount < cyclicCount:
                        group.NAME = 'Cyclic' + group.NAME
        return 

    def determineAlcoholGroups(self, functionalGroups):
        """ Appends the functional groups created from alocholic oxygens, with ring classifications, to the functionalGroups list passed in.

            functionalGroups (list) : List of Molecule obejects that represent the functional groups of a SMILES code

            Notes:
                This function looks to see if the alcoholic oxygen is attached within an aromatic, non-aromatic ring strucutres, or no at all
                Build upon alochol indices found in the molecule
        """

        for index in self.ALCOHOLICINDICES:

            alohcolBond = self.bondData[index]
            alochol = Molecule("RO", "Alcohol")
            alochol.atomData[0].index = alohcolBond[0].index
            alochol.atomData[1].index = index

            if alohcolBond[0].index in self.AROMATICINDICES:
                alochol.NAME = "AromaticAlcohol"
            elif alohcolBond[0].index in self.CYCLICINDICES:
                alochol.NAME = "CyclicAlcohol"

            functionalGroups.append(alochol)

        return 

    def priorityFilter(self, functionalGroups):
        """ Filters functionalGroups list of all non-maximized R groups representations of functional group families

            functionalGroups (list) : List of Molecule obejects that represent the functional groups in a smiles code

            For exmaple, a secondary amine and a tertiary amine. The tertiary amine is the more complete strucuture.
            
            RN(R)R vs RNR

            Heirarchy scrub will remove the lower of the two from these relations, but keep non-R group contained ones such as ketone and ester.

            RC(=O) vs RC(=O)O

            Because the number of R groups is the same, this group will not be filtered, and instead by considered as its own individual group
        """

        index = -1
        while index != len(functionalGroups) - 1:

            index += 1
            group = functionalGroups[index]
            groupAtoms = group.atomData
            mainGroupAtoms = self.filterRgroup(groupAtoms)
            numMainRAtoms = len(self.getRgroups(groupAtoms))

            for compareIndex, compareGroup in enumerate(functionalGroups):

                compareAtoms = compareGroup.atomData
                mainCompareAtoms = self.filterRgroup(compareAtoms)
                numCompareRAtoms = len(self.getRgroups(compareAtoms))

                if(
                    all(i in mainCompareAtoms for i in mainGroupAtoms) 
                    and all(i in mainGroupAtoms for i in mainCompareAtoms) 
                    and compareIndex != index
                ):
                    if numMainRAtoms > numCompareRAtoms:
                        del(functionalGroups[compareIndex])
                        index = -1
                        break
                    elif numCompareRAtoms > numMainRAtoms:
                        del(functionalGroups[index])
                        index = -1
                        break
        return 

    def findPreciseGroups(self, functionalGroups):
        """ Filters functional groups for their precise groups.

            functionalGroups (list) : List of Molecule obejects that represent the functional groups in a smiles cod.

            Relations such as ketone in ester, amine in amide, ether in ester are elimated. 
        """
        preciseFgs = functionalGroups[:]
        index = -1
        while index != len(preciseFgs) - 1:

            index += 1
            group = preciseFgs[index]
            groupAtoms = group.atomData
            groupIndices = []
            for atom in groupAtoms.values():
                groupIndices.append(atom.index)
            for compareIndex, compareGroup in enumerate(preciseFgs):

                if compareIndex == index:
                    continue

                compareAtoms = compareGroup.atomData
                compareIndices = []
                for atom in compareAtoms.values():
                    compareIndices.append(atom.index)

                if len(compareIndices) != len(groupIndices):

                    if all(i in groupIndices for i in compareIndices):
                        del(preciseFgs[compareIndex])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                    elif all(i in compareIndices for i in groupIndices):
                        del(preciseFgs[index])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                elif all(i in groupIndices for i in compareIndices) and all(i in compareIndices for i in groupIndices):

                    numMainRAtoms = len(self.getRgroups(groupAtoms))
                    numCompareRAtoms = len(self.getRgroups(compareAtoms))

                    if numMainRAtoms > numCompareRAtoms:
                        del(preciseFgs[index])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                    elif numCompareRAtoms > numMainRAtoms:
                        del(preciseFgs[compareIndex])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                del(compareIndices)

        return preciseFgs

    def detetminePrimaryAmine(self, nitrogenIndex):
        """ Return primary amine Molecule object with proper ring classification, if it exists on a ring
            Return false if nitrogen is not a primary amine

            nitrogenIndex (int) : The index of the nitrogen in the smiles code to be analyzed for being a primary amine

        """

        # Check bonds on nitrogen to see if it is a Primary Amine
        smilesBonds = self.bondData[nitrogenIndex]

        # Primary amines have one, single bonded partner (i.e. no explicit bonds =/#), which may be in a ring
        if len(smilesBonds) == 1 and not self.BOND_REGEX.findall(smilesBonds[0].symbol):

            bondedIndex = smilesBonds[0].index

            primaryAmine = Molecule("RN", "PrimaryAmine")

            # Set the R group index to that of the smiles bonded nitrogen index
            primaryAmine.atomData[0].index = bondedIndex
            # Set the nitrogen index to that of the smiles code nitrogen index
            primaryAmine.atomData[1].index = nitrogenIndex

            if bondedIndex in self.AROMATICINDICES:
                primaryAmine.NAME = "AromaticPrimaryAmine"

            elif bondedIndex in self.CYCLICINDICES:
                primaryAmine.NAME = "CyclicPrimaryAmine"

            return primaryAmine

        return False
