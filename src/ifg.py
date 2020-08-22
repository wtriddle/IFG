from Molecule import Molecule
import re


class ifg(Molecule):

    def __init__(self, smiles, REFCODE):
        super().__init__(smiles, REFCODE)
        self.functionalGroups = self.findFunctionalGroups()
        self.preciseFunctionalGroups = self.findPreciseGroups(
            self.functionalGroups)

    def findFunctionalGroups(self):

        functionalGroups = []

        for atom in self.atomData:
            groups = self.whichGroup(atom)

            if atom[1] == 'N':
                primaryAmine = self.detetminePrimaryAmine(atom)
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
        symbol = atom[1]
        index = atom[0]

        for line in open('FGlist.txt', 'r'):
            lineInfo = re.compile(r'\S+').findall(line)
            lineInfo[0] = lineInfo[0].replace('[R]', 'R')
            template = Molecule(lineInfo[0], lineInfo[1])

            if len(self.chargeRegex.findall(symbol)) == 0 and len(self.chargeRegex.findall(template.SMILES)) != 0:
                continue

            if template.ALCOHOLICINDICES and index not in self.ALCOHOLICINDICES:
                continue

            if symbol in template.SMILES:
                expansionPoint = 0

                # R groups from temple.bondData can be removed because their expansion leads back to a main group atom from which is was expnaded from
                for tempAtom in template.atomData:
                    if tempAtom[1] == 'R':
                        template.bondData[tempAtom[0]].clear()

                # Find expansion indicies or points inside template
                for tempAtom in template.atomData:

                    if not template.ALCOHOLICINDICES:
                        if symbol == tempAtom[1]:
                            expansionPoint = tempAtom[0]
                            break
                    else:
                        expansionPoint = template.ALCOHOLICINDICES[0]
                        break

                # Set up for expansion call at expansionPoint
                template.atomData[expansionPoint][0] = atom[0]
                self.removeBond(expansionPoint, template)

                expandGroup = self.expandGroup(
                    atom, expansionPoint, template, None, atom[0])
                if expandGroup:
                    matches.append(template)

        return matches

    def expandGroup(self, atom, expansionPoint, template, skipIndex=None, smilesIndexInit=None):

        smilesIndices = []  # Smiles indicies which have already been used in template
        smilesIndices.append(smilesIndexInit)
        templateBonds = template.bondData[expansionPoint]
        mainTemplateBonds = self.getMainGroups(templateBonds)
        smilesBonds = self.bondData[atom[0]]
        templateAtoms = template.atomData

        # Not enough bonds at desposal of smiles code to satisfy template means False match
        if len(templateBonds) > len(smilesBonds):
            return False

        # Main group Template Analysis
        for tempBond in mainTemplateBonds:
            tempIndex = tempBond[0]
            if tempIndex == skipIndex:
                continue
            tempSymbol = tempBond[1]
            for smilesBond in smilesBonds:

                smilesIndex = smilesBond[0]
                smilesSymbol = smilesBond[1]

                if tempSymbol == smilesSymbol and smilesIndex not in smilesIndices:

                    path = self.expandGroup(
                        smilesBond, tempIndex, template, expansionPoint, atom[0])

                    if path:
                        smilesIndices.append(smilesIndex)
                        templateAtoms[tempIndex][0] = smilesIndex
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

                    smilesIndex = smilesBond[0]
                    smilesSymbol = smilesBond[1]

                    if smilesIndex in smilesIndices or smilesSymbol[0] in self.BONDS:
                        continue
                    RgroupCounter += 1
                    if RgroupCounter < len(Rgroups) - 1:
                        templateAtoms[Rgroups[RgroupCounter]
                                      [0]][0] = smilesIndex
                    else:
                        templateAtoms[Rgroups[RgroupCounter]
                                      [0]][0] = smilesIndex
                        break
                else:
                    return False

                return True

            else:
                return False

    def getRgroups(self, atomSet):

        Rgroups = []

        for atom in atomSet:
            if atom[1] == 'R':
                Rgroups.append(atom)

        return Rgroups

    def getMainGroups(self, atomSet):

        mainGroups = []

        for atom in atomSet:
            if atom[1] != 'R':
                mainGroups.append(atom)

        return mainGroups

    def removeBond(self, index, template):

        for i, bonds in enumerate(template.bondData):
            for j, bond in enumerate(bonds):
                if bond[0] == index:
                    del(template.bondData[i][j])
        return 0

    def repetitionScrub(self, functionalGroups):

        index = -1

        while index != len(functionalGroups) - 1:

            index += 1
            group = functionalGroups[index]
            groupAtoms = group.atomData

            for compareIndex, compareGroup in enumerate(functionalGroups):

                compareAtoms = compareGroup.atomData

                if all(i in compareAtoms for i in groupAtoms) and all(i in groupAtoms for i in compareAtoms) and compareIndex != index:
                    del(functionalGroups[compareIndex])
                    index = - 1
                    break

        return 0

    def determineCyclicGroups(self, functionalGroups):

        for group in functionalGroups:
            groupAtoms = group.atomData
            index = aromaticCount = cyclicCount = 0
            while index < len(groupAtoms) - 1:
                indicies = [groupAtoms[index][0], groupAtoms[index+1][0]]
                if all(i in self.AROMATICINDICES for i in indicies):
                    aromaticCount += 1
                elif all(i in self.CYCLICINDICES for i in indicies):
                    cyclicCount += 1
                index += 1
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
            alochol = molecule("RO", "Alcohol")
            alochol.atomData[0][0] = alohcolBond[0][0]
            alochol.atomData[1][0] = index

            if alohcolBond[0][0] in self.AROMATICINDICES:
                alochol.NAME = "AromaticAlcohol"
            elif alohcolBond[0][0] in self.CYCLICINDICES:
                alochol.NAME = "CyclicAlcohol"

            functionalGroups.append(alochol)

        return 0

    def heirarchyScrub(self, functionalGroups):

        index = -1
        while index != len(functionalGroups) - 1:

            index += 1
            group = functionalGroups[index]
            groupAtoms = group.atomData
            mainGroupAtoms = self.getMainGroups(groupAtoms)
            numMainRAtoms = len(self.getRgroups(groupAtoms))

            for compareIndex, compareGroup in enumerate(functionalGroups):

                compareAtoms = compareGroup.atomData
                mainCompareAtoms = self.getMainGroups(compareAtoms)
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
            for atom in groupAtoms:
                groupIndices.append(atom[0])
            for compareIndex, compareGroup in enumerate(preciseFgs):

                if compareIndex == index:
                    continue

                compareAtoms = compareGroup.atomData
                compareIndices = []
                for atom in compareAtoms:
                    compareIndices.append(atom[0])

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

    def detetminePrimaryAmine(self, atom):

        smilesBonds = self.bondData[atom[0]]

        if len(smilesBonds) == 1 and smilesBonds[0][1][0] not in self.BONDS:

            primaryAmine = molecule("RN", "PrimaryAmine")
            primaryAmine.atomData[0][0] = smilesBonds[0][0]
            primaryAmine.atomData[1][0] = atom[0]
            if smilesBonds[0][0] in self.AROMATICINDICES:
                primaryAmine.NAME = "AromaticPrimaryAmine"
            elif smilesBonds[0][0] in self.CYCLICINDICES:
                primaryAmine.NAME = "CyclicPrimaryAmine"
            return primaryAmine

        else:
            return False
