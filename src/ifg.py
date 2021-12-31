""" Algorithm which analyzes the atom and bond data of a given SMILES code to produce counts of pre-defined organic functional groups in FGlist.txt.

    This algorithm is built from the Molecule class, which has two important data structures: atomData and bondData. These two
    dictionaries possess the data of how the molecule atoms are connected. bondData contains all the possible bond paths that stem from
    every atom in the SMILES code. The same is true of decoded functional group templates from the Molecule class.
    Therefore, if a functional group template can model itself upon at least one of those paths in the SMILES code, 
    then that functional group must be present inside of the SMILES code. The concept of a valid path for a functional group
    within the SMILES code string is the solution to this algorithm.

Key Attributes:

    The following two functional group data have keys for the counts of their that specific SMILES code functional groups, 
    including extra classification like cyclic or aromatic. Also includes ring data and alcohol counts in both. They differ slightly:

    all_fgs (Dict) : The set of FGS after.
    Includes overlapping functional groups, such as a ketone inside of an ester. 
    
    exact_fgs (Dict) : The dictionary which consider overlapping funcitonal groups.
    Excludes overlapping functional groups, ex ketones left unrecorded when seen inside esters, etc.
"""


from Molecule import Molecule
from constants import BONDS, CHARGE_REGEX, NON_BRANCHING_SYMBOLS, BOND_REGEX
from helpers import formatSmiles
import re
from config import FGSPATH
from datetime import datetime


class ifg(Molecule):
    """ The class algorithm which identifies the functional group counts of a SMILES code. """

    def __init__(self, SMILES, REFCODE):
        """ Determine the functional groups from a SMILES code:
            
            ## Parameters:
            SMILES (string) : A valid simplified molecular input line entry system (SMILES) code
            REFCODE (string) : The referenced code for this particular SMILES code

            ## Usage:
            >>> SMILES = "OC(=O)CC1CCC(=O)CC1"
            >>> REFCODE = "KUZQIG"
            >>> fgs = ifg(SMILES, REFCODE)
            >>> print(fgs)
            all_fgs: [CarboxylicAcid, Ketone, CyclicKetone, Alcohol]
            exact_fgs: [CarboxylicAcid, CyclicKetone]
            >>> data_all = createDataDict(fgs.all_fgs)
            >>> data_precise = createDataDict(fgs.exact_fgs)
            >>> ...
        """

        super().__init__(SMILES, REFCODE)                           # Create Molecule object from input SMILES code
        self.dfs_time = datetime.now()                              # Performance metric for DFS portion of the algorithm (main IFG)
        self.filter_time = datetime.now()                           # Performance metric for the filter portion of the algorithm
        self.all_fgs = self.findFunctionalGroups()                  # all_fgs is a property ready after IFG
        self.exact_fgs = self.overlapFilter(self.all_fgs)           # Filter overlapping groups to get different format of the data

    def findFunctionalGroups(self):
        """ Return a list of Molecule objects which represent the functional groups of a particular SMILES code

            Notes:
                FGs as Molecules inside of the SMILES allows the atomic indicies to
                reflect the positioning within the SMILES
            """

        # List of identified functional groups
        matches_list = []

        now = datetime.now()
        # Main IFG loop (loops for num_atoms*num_fgs times)
        for atom in self.atomData.values():                     # Loop over all atoms in the SMILES
            groups = self.whichGroup(atom)                      # Determine the functional groups which stem from an atom (loops over all FGS)

            if atom.symbol == 'N':                              # Special handling for nitrogens
                primaryAmine = self.detetminePrimaryAmine(atom.index)

                if primaryAmine:
                    matches_list.append(primaryAmine)

            for group in groups:                                # For each group (each are Molecule objects) found to stem from an atom
                matches_list.append(group)                      # Add it to the "Matches list" of identified FG
        
        self.dfs_time = datetime.now() - now

        # Post data processing and classification methods
        now = datetime.now()
        self.repetitionFilter(matches_list)                 # Remove repeated groups
        self.hierarchyFilter(matches_list)                  # Remove hierarchically overlapped groups

        self.classifyCyclic(matches_list)            # Identify the groups which are aromatic and non-aromatic
        self.determineAlcoholGroups(matches_list)           # Create Molecule objects for identified alcohols
        self.filter_time = datetime.now() - now
        return matches_list

    def whichGroup(self, atom):
        """ Returns the list of functional groups as Molecule objects that stem from a given atom in the SMILES code

            atom (type Atom) : An atom object who is to be analyzed for being contained in a possible functional group

            Process of whichGroup:
                Filters functional groups which do not contain given atom symbol
                Selects an expansion point from within a potential functional group to start expansion from 
                    (i.e. Which C to choose in [R]C(=O)NC(=O)[R], 1st or 2nd)
                Creates Molecule object of functional group
                Calls expandgroup to complete path validation 
            
            Notes:
                Loops over every possible functional group in FGlist
                Every atom has a possibility to be branched into a functional group. Therefore, every atom must be looked at individually.
                R group paths are elimated because they cannot be expanded from to validate the functional group. Only main group atoms can be expanded from
        """

        matches = []
        for (FGsmiles, name) in [                                                   # Fetch the (smiles, name) functional group pairs from FGlist.txt
            [y.strip() for y in x.split(' ')]                                       # Split into a formatted (smiles, name) pair
            for x                                                                   # Done for each line
            in open(FGSPATH.resolve(), "r+").readlines()                            # Create list using all lines in FGlist.txt
        ]:
            FGsmiles = formatSmiles(FGsmiles).replace('[R]', 'R')                   # Remove [R] from brackets

            if(                                                                     # Non-charged SMILES cannot contain a charged FG
                len(CHARGE_REGEX.findall(atom.symbol)) == 0                         # Charged FG's
                and len(CHARGE_REGEX.findall(FGsmiles)) != 0                        # In a non-charged SMILES
            ):
                continue                                                            # Are skipped

            if atom.symbol in FGsmiles:                                             # If the atom symbol argument is inside the FGtemplate
                FGtemplate = Molecule(FGsmiles, name)                               # FG Molecule object

                if(                                                                 # Alcoholic FG's can only stem from alcoholic oxygens
                    FGtemplate.ALCOHOLICINDICES                                     # Alcoholic FG
                    and atom.index not in self.ALCOHOLICINDICES                     # With a Non-alcoholic oxygen atom
                ):
                    continue                                                        # Are skipped

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
                        expansionPoint = FGtemplate.ALCOHOLICINDICES[0]             # Choose 0th one (if FG exists in molecule, match will occur)
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

    def expandGroup(self, atom, expansionPoint, template, skipIndex=None, validSmilesIndex=None):
        """ Return true if an atom can be expanded into the given template

            atom (type Atom) : the SMILES code atom which is to be analyzed for having a valid bond with respect to the template
            expansionPoint (int) : The position inside of the template for which atom tested as having the matching bond paths inside of the template with respect to the SMILES code bonding paths
            template (type Molecule) : The functional group for which the given atom will be checked against.
            skipIndex (int, default=None) : A template bond for which to skip over
            validSmilesIndex (int, default=None) : The indices confirmed in the SMILES code to be apart of a functional group, with designated positions
        """

        smilesIndices = []                                      # Smiles indicies which have already been mapped into the template 
        smilesIndices.append(validSmilesIndex)                  # Add the validated SMILES atom index to the list of mapped SMILES to FG indicies
        templateBonds = template.bondData[expansionPoint]       # Atoms which stem from the atom of interest within FGtemplate
        mainTemplateBonds = self.filterRgroup(templateBonds)    # Main group atoms which stem from atom of interest within FGtemplate
        smilesBonds = self.bondData[atom.index]                 # Atoms which stem from the atom of interst within the SMILES 
        templateAtoms = template.atomData                       # All atoms inside of the FGtemplate

        if len(templateBonds) > len(smilesBonds):               # Not enough bonds from SMILES atom to satisfy FGtemplate atom means False match
            return False

        for tempBond in mainTemplateBonds:                      # Main group atom analysis of FGtemplate Analysis
            tempIndex = tempBond.index                          
            if tempIndex == skipIndex:                          # Skip the index which was just stemmed from
                continue
            tempSymbol = tempBond.symbol               

            for smilesBond in smilesBonds:                      # For all bonding paths in the SMILES

                smilesIndex = smilesBond.index                 
                smilesSymbol = smilesBond.symbol

                if(
                    tempSymbol == smilesSymbol                  # Same bonding symbol in SMILES and FGTemplate is found
                    and smilesIndex not in smilesIndices        # And the SMILES atom bond is unsed in FG bonding path
                ):                                              # Means a valid bonding path has been located in the SMILES with respect to the FGtemplate

                    path = self.expandGroup(smilesBond,         # Use SMILES bonded atom as new point to expand from
                                            tempIndex,          # Point of expansion in template is the atom index in template which matched
                                            template,           # Use same template object for all recurisve calls
                                            expansionPoint,     # Expansion point of previous path is now skipped since it was an atom that was just validated on this function call
                                            atom.index          # SMILES atom index is added to list of validated atoms
                    )                                           # New expansion path is now set with these arguments

                    if path:                                    # If a bonding path is valid
                        smilesIndices.append(smilesIndex)       # Add the current SMILES index to list of valid atoms
                        templateAtoms[tempIndex].index = smilesIndex    # Set the index of the FGtemplate to reflect its mapped SMILES index
                        break
                    else:                                       # No path in most recent direction means try another path
                        continue                                
            else:                                               # No SMILES bonds to fulfill path means invalid FG in SMILES
                return False    

        else:                                                   # Rgroup FGtemplate Analysis

            RgroupCounter = -1                                  # Loop counter for Rgroups
            Rgroups = self.getRgroups(templateBonds)            # Filter template bonds for their Rgropus

            if not Rgroups:                                     # No Rgroups after main group analysis is a valid FG
                return True

            if len(smilesBonds) > len(mainTemplateBonds):       # SMILES bonds must have at least 1 more bond than the main group atom to fulfill Rgroup requirement

                for smilesBond in smilesBonds:                  # Loop over all SMILES bonds which stem from atom

                    smilesIndex = smilesBond.index
                    smilesSymbol = smilesBond.symbol

                    if(
                        smilesIndex in smilesIndices            # A used SMILES atom index
                        or smilesSymbol[0] in BONDS        # Or a non single-bonded atom
                    ):                                          # Are invalid for R group requirements, skip to next bonds
                        continue


                    RgroupCounter += 1                          # Go to next Rgroup
                    if RgroupCounter < len(Rgroups) - 1:        # If there are Rgroups to process
                        templateAtoms[                          # Set index of Rgroup atom in template
                            Rgroups[RgroupCounter].index        # Identified by the Rgroup list
                        ].index = smilesIndex                   # As the SMILES code atom singly bonded to the atom of interest
                    else:                                       # Process final R groups and break after setting index
                        templateAtoms[
                            Rgroups[RgroupCounter].index
                        ].index = smilesIndex
                        break
                else:                                           # No SMILES bonds available means invalid FG
                    return False

                return True                                     # If all main group atoms and Rgroup atom requirements are met, then a valid FG has been located

            else:                                               # Not enough atoms in SMILES to fulfill R group requirements
                return False

    def getRgroups(self, atomSet):
        """ Returns a dictionary or list of all R group atoms in a given atomSet.

            atomSet (list || dict) : Data container which has Atom objects, either as values or in a list

            Notes:
                A dictionary is keyed by index and valued by its atom symbol
        """

        if isinstance(atomSet, list):
            return [atom for atom in atomSet if atom.symbol == 'R']
        if isinstance(atomSet, dict):
            return {atom.index: atom for atom in atomSet.values() if atom.symbol == 'R'}

    def filterRgroup(self, atomSet):
        """ Returns a dictionary or list of all non R group atoms in a given atomSet.

            atomSet (list || dict) : Data container which has Atom objects, either as values or in a list

            Notes:
                A dictionary is keyed by index and valued by its atom symbol
        """

        if isinstance(atomSet, list):
            return [atom for atom in atomSet if atom.symbol != 'R']
        if isinstance(atomSet, dict):
            return {atom.index: atom for atom in atomSet.values() if atom.symbol != 'R'}

    def removeBond(self, index, template):
        """ Removes a bond path from template bondData

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

    def classifyCyclic(self, functionalGroups):
        """ Label functional groups with ring classification based on the ring structure of Molecule (self) 

            functionalGroups (list) : List of Molecule obejects that represent the functional groups in a smiles code

            Notes: 
                Some rings are directly connected to their opposite type, i.e aromatic-nonAromatic
                To be the most accurate, a tally system was created to tally how many valid bonds fall under a certain classification.
                Whichever ring type has more tallies, that one is taken as the correct classification
        """

        for group in functionalGroups:                              # Loop over all FG's, group is a Molecule object
            groupAtoms = group.atomData                             # Atoms part of the FG with SMILES validated indicies
            aromaticCount = cyclicCount = 0                         # Tallies for atoms part of a specific ring
            for templateIndex, smilesAtom in groupAtoms.items():    # Loop over all atoms, with their indicies, in the FG

                if smilesAtom.symbol not in NON_BRANCHING_SYMBOLS:     # Non-linear symbols are never cyclic
                    continue

                if(                                                 
                    smilesAtom.index not in self.AROMATICINDICES    # Non-aromatic
                    and smilesAtom.index not in self.CYCLICINDICES  # And non-cyclic atoms
                ):
                    continue                                        # Cannot be part of any rings, so skip these atoms

                for templateBond in group.bondData[templateIndex]:          # Loop over all bonds which stem from a given FG atom
                    smilesBond = groupAtoms[templateBond.index]             # Find the SMILES bond equivalent of the FG bond
                    indicies = [smilesAtom.index, smilesBond.index]         # Concatinate the atom and bonded atom indicies together
                    if all(i in self.AROMATICINDICES for i in indicies):    # If both atoms are part of an aromatic ring
                        aromaticCount += 1                                  # Add to aromatic tally
                    elif all(i in self.CYCLICINDICES for i in indicies):    # If both atoms are part of a cyclic ring
                        cyclicCount += 1                                    # Add to cyclic tally
            else:                                                   # Once the FG has been evaluated for aromatic/cyclic tallies
                if aromaticCount or cyclicCount:                    # Validate that there were aromatic/cyclic tallies
                    if aromaticCount >= cyclicCount:                # Take aromatic nomenclature if more tallies for aromatic. Equality is deffered to aromatic
                        group.NAME = 'Aromatic' + group.NAME
                    elif aromaticCount < cyclicCount:               # Otherwise use cyclic nomenclature if more tallies for cyclic
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

    def hierarchyFilter(self, functionalGroups):
        """ Filters functionalGroups list of all non-maximized R groups representations of functional group families

            functionalGroups (list) : List of Molecule objects that represents the functional groups in a smiles code

            Notes:
                For exmaple, a secondary amine and a tertiary amine. The tertiary amine is the more complete strucuture.
                RN(R)R vs RNR
                Uses a comparsion for R group with more R atoms given that the non-R group structures are equal 
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

    def overlapFilter(self, functionalGroups):
        """ Filters functional groups list based on full overlapping condition

            functionalGroups (list) : List of Molecule obejects that represent the functional groups in a SMILES

            Notes:
                Relations such as ketone in ester, amine in amide, ether in ester are elimated using this function
        """
        exact_fgs = functionalGroups[:]
        index = -1
        while index != len(exact_fgs) - 1:

            index += 1
            group = exact_fgs[index]
            groupAtoms = group.atomData
            groupIndices = []
            for atom in groupAtoms.values():
                groupIndices.append(atom.index)
            for compareIndex, compareGroup in enumerate(exact_fgs):

                if compareIndex == index:
                    continue

                compareAtoms = compareGroup.atomData
                compareIndices = []
                for atom in compareAtoms.values():
                    compareIndices.append(atom.index)

                if len(compareIndices) != len(groupIndices):

                    if all(i in groupIndices for i in compareIndices):
                        del(exact_fgs[compareIndex])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                    elif all(i in compareIndices for i in groupIndices):
                        del(exact_fgs[index])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                elif(
                    all(i in groupIndices for i in compareIndices) 
                    and all(i in compareIndices for i in groupIndices)
                ):

                    numMainRAtoms = len(self.getRgroups(groupAtoms))
                    numCompareRAtoms = len(self.getRgroups(compareAtoms))

                    if numMainRAtoms > numCompareRAtoms:
                        del(exact_fgs[index])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                    elif numCompareRAtoms > numMainRAtoms:
                        del(exact_fgs[compareIndex])
                        del(groupIndices)
                        del(compareIndices)
                        index = -1
                        break

                del(compareIndices)

        return exact_fgs

    def detetminePrimaryAmine(self, nitrogenIndex):
        """ Return primary amine Molecule object with proper ring classification, if it exists on a ring
            Return false if nitrogen is not a primary amine

            nitrogenIndex (int) : The index of the nitrogen in the smiles code to be analyzed for being a primary amine

        """

        nitrogenBonds = self.bondData[nitrogenIndex]                    # Get bonds connected to nitrogen

        if(
            len(nitrogenBonds) == 1                                     # Primary amines have one bond
            and not BOND_REGEX.findall(nitrogenBonds[0].symbol)    # That must be a single bond
        ):                                                              # Otherwise, nitrogen is not a primary amine

            bondedIndex = nitrogenBonds[0].index                        

            primaryAmine = Molecule("RN", "PrimaryAmine")               # Create primary amine Molecule object

            primaryAmine.atomData[0].index = bondedIndex                # Set the R group index to the smiles bonded nitrogen index
            primaryAmine.atomData[1].index = nitrogenIndex              # Set the nitrogen index to the smiles code nitrogen index

            if bondedIndex in self.AROMATICINDICES:                     # Add aromatic nomenclature if necessary
                primaryAmine.NAME = "AromaticPrimaryAmine"

            elif bondedIndex in self.CYCLICINDICES:                     # Or cyclic nomenclature 
                primaryAmine.NAME = "CyclicPrimaryAmine"

            return primaryAmine                                         # Return the primary amine object if true

        return False                                                    # Otherwise return false

    def printall_fgs(self):
        for f in self.all_fgs:
            print(f)
    def printSpecificFgs(self, fgs):
        for f in fgs:
            print(f, f.getSymbolDict())

    def __str__(self):
        """ String function for printing out fgs object to screen
            >>> fgs = ifg(SMILES, REFCODE)
            >>> print(fgs)
            __str__ output
        """

        all_fgs = self.all_fgs
        exact_fgs = self.exact_fgs

        # All FGS list to screen
        s = "all_fgs: ["
        for f in all_fgs:
            s+=f.NAME
            s+=", "
        s=s[0:-2]       # , correction
        s+=']\n'

        # EXACT FGS list to screen
        s+= "exact_fgs: ["
        for f in exact_fgs:
            s+=f.NAME
            s+=", "
        s=s[0:-2]       # , correction
        s+=']\n'

        # Display constructed string
        return s