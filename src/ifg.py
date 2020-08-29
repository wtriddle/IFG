from Molecule import Molecule
import re
from helpers import createFgDataDict
import os


class ifg(Molecule):

    def __init__(self, smiles, REFCODE):
        super().__init__(smiles, REFCODE)
        self.allFgs = self.findFunctionalGroups()
        self.preciseFgs = self.findPreciseGroups(
            self.allFgs)

    def findFunctionalGroups(self):
        """ Return a list of Molecule objects, where the indices of the atoms in the Molecule reflect those inside of the smiles code

        """
        functionalGroups = []

        for atom in self.atomData.values():
            groups = self.whichGroup(atom)

            if atom.symbol == 'N':
                # Pass in nitrogen index to see if nitrogen is a Primary Amine
                primaryAmine = self.detetminePrimaryAmine(atom.index)
                if primaryAmine:
                    functionalGroups.append(primaryAmine)

            for group in groups:
                functionalGroups.append(group)

        self.repetitionScrub(functionalGroups)
        self.heirarchyScrub(functionalGroups)
        self.determineCyclicGroups(functionalGroups)
        self.determineAlcoholGroups(functionalGroups)

        return functionalGroups

    def whichGroup(self, atom):
        matches = []
        for line in open(os.getcwd() + '/src/resources/FGlist.txt', 'r'):
            lineInfo = re.compile(r'\S+').findall(line)
            lineInfo[0] = lineInfo[0].replace('[R]', 'R')
            template = Molecule(lineInfo[0], lineInfo[1])

            # Skip charged FG's unless the atom is charged
            if len(self.CHARGE_REGEX.findall(atom.symbol)) == 0 and len(self.CHARGE_REGEX.findall(template.SMILES)) != 0:
                continue

            # Skip alcohol-containing FG's unless the atom is an alcohol (remember the H is implicit)
            if template.ALCOHOLICINDICES and atom.index not in self.ALCOHOLICINDICES:
                continue

            # If the symbol is inside the FG template
            if atom.symbol in template.SMILES:
                expansionPoint = 0

                # R groups from temple.bondData can be removed because their expansion leads back to a main group atom from which is was expnaded from
                for tempAtom in template.atomData.values():
                    if tempAtom.symbol == 'R':
                        template.bondData[tempAtom.index].clear()

                # Find expansion indicies or points inside template
                for tempAtom in template.atomData.values():

                    if not template.ALCOHOLICINDICES:
                        if atom.symbol == tempAtom.symbol:
                            expansionPoint = tempAtom.index
                            break
                    else:
                        expansionPoint = template.ALCOHOLICINDICES[0]
                        break

                # Set up for expansion call at expansionPoint
                template.atomData[expansionPoint].index = atom.index
                self.removeBond(expansionPoint, template)

                expandGroup = self.expandGroup(
                    atom, expansionPoint, template, None, atom.index)
                if expandGroup:
                    matches.append(template)

        return matches

    def expandGroup(self, atom, expansionPoint, template, skipIndex=None, smilesIndexInit=None):
        """
        atom.symbol
        atom.index
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
                        #
                        templateAtoms[Rgroups[RgroupCounter].index
                                      ].index = smilesIndex
                    else:
                        #
                        templateAtoms[Rgroups[RgroupCounter].index
                                      ].index = smilesIndex
                        break
                else:
                    return False

                return True

            else:
                return False

    def getRgroups(self, atomSet):
        if isinstance(atomSet, list):
            return [atom for atom in atomSet if atom.symbol == 'R']
        if isinstance(atomSet, dict):
            return {atom.index: atom for atom in atomSet.values() if atom.symbol == 'R'}

    def filterRgroup(self, atomSet):
        if isinstance(atomSet, list):
            return [atom for atom in atomSet if atom.symbol != 'R']
        if isinstance(atomSet, dict):
            return {atom.index: atom for atom in atomSet.values() if atom.symbol != 'R'}

    def removeBond(self, index, template):

        for i, bonds in enumerate(template.bondData.values()):
            for j, bond in enumerate(bonds):
                if bond.index == index:
                    del(template.bondData[i][j])
        return 0

    def repetitionScrub(self, functionalGroups):

        index = -1

        while index != len(functionalGroups) - 1:

            index += 1
            group = functionalGroups[index]
            groupAtoms = group.atomData.values()

            for compareIndex, compareGroup in enumerate(functionalGroups):

                compareAtoms = compareGroup.atomData.values()

                if all(i in compareAtoms for i in groupAtoms) and all(i in groupAtoms for i in compareAtoms) and compareIndex != index:
                    del(functionalGroups[compareIndex])
                    index = - 1
                    break

        return 0

    def determineCyclicGroups(self, functionalGroups):

        for group in functionalGroups:
            groupAtoms = group.atomData
            # Only start from main group atoms for cyclic classification
            aromaticCount = cyclicCount = 0
            for templateIndex, smilesAtom in groupAtoms.items():

                if smilesAtom.symbol not in self.LINEARSYMBOLS:
                    continue

                # If maingroup atom does not have any cyclic property, then skip it
                if smilesAtom.index not in self.AROMATICINDICES and smilesAtom.index not in self.CYCLICINDICES:
                    continue

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
        return 0

    def determineAlcoholGroups(self, functionalGroups):

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

        return 0

    def heirarchyScrub(self, functionalGroups):

        index = -1
        while index != len(functionalGroups) - 1:

            index += 1
            group = functionalGroups[index]
            groupAtoms = group.atomData
            mainGroupAtoms = self.filterRgroup(groupAtoms)
            # print(groupAtoms)
            numMainRAtoms = len(self.getRgroups(groupAtoms))

            for compareIndex, compareGroup in enumerate(functionalGroups):

                compareAtoms = compareGroup.atomData
                mainCompareAtoms = self.filterRgroup(compareAtoms)
                numCompareRAtoms = len(self.getRgroups(compareAtoms))

                if all(i in mainCompareAtoms for i in mainGroupAtoms) and all(i in mainGroupAtoms for i in mainCompareAtoms) and compareIndex != index:
                    if numMainRAtoms > numCompareRAtoms:
                        del(functionalGroups[compareIndex])
                        index = -1
                        break
                    elif numCompareRAtoms > numMainRAtoms:
                        del(functionalGroups[index])
                        index = -1
                        break
        return 0

    def findPreciseGroups(self, functionalGroups):

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
