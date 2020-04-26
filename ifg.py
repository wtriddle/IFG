from molecule import molecule
import re

class ifg(molecule):

	def __init__(self,smiles, REFCODE):
		super().__init__(smiles, REFCODE)
		self.functionalGroups = self.findFunctionalGroups()
	
	def findFunctionalGroups(self):
	
		functionalGroups = []

		for atom in self.atomData:
			groups = self.whichGroup(atom)

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
		# if index in self.ALCOHOLICINDICES:
		# 	print("analyzing alcoholic", atom)

		for line in open('FGlist.txt','r'):
			lineInfo = re.compile(r'\S+').findall(line)
			lineInfo[0] = lineInfo[0].replace('[R]','R')
			template = molecule(lineInfo[0], lineInfo[1])

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

				expandGroup = self.expandGroup(atom,expansionPoint,template, [])
				if expandGroup:
					if template.SMILES == 'RC(=O)N(R)C(=O)R':
						print("\nMatched ", template.SMILES)
						print(template.atomData)
						print('\n')
					matches.append(template)
				else:
					if template.SMILES == 'RC(=O)N(R)C(=O)R':
						print("expansion Failse")
		
		return matches

	def expandGroup(self, atom, expansionPoint, template, smilesIndexMatches):
		
		atomMatches = []
		RgroupMatches = []
		tempIndexMatches = []

		templateBonds = template.bondData[expansionPoint]
		templateAtoms = template.atomData
		Rgroups = self.getRgroups(templateBonds)
		RgroupCount = -1
		smilesBonds = self.bondData[atom[0]]
		if len(templateBonds) > len(smilesBonds):
			if template.SMILES == 'RC(=O)N(R)C(=O)R':
				print("templateBonds = ", templateBonds)
				print("smilesBonds = ", smilesBonds)
			return False
		if template.SMILES == 'RC(=O)N(R)C(=O)R':
			print("\nInitial Conditions:")
			print("templateSmiles = ", template.SMILES)
			print("smilesAtom = ", self.atomData[atom[0]])
			print("smilesBonds = ",  smilesBonds)
			print("expansionPoint = ", expansionPoint)
			print("templateBonds = ", templateBonds)
			print("templateAtoms = ", templateAtoms)
			print("Rgroups = ", Rgroups)

		templateAtoms[expansionPoint][0] = atom[0]
		smilesIndexMatches.append(atom[0])
		self.removeBond(expansionPoint,template)

		for smilesBond in smilesBonds:

			index = smilesBond[0]
			symbol = smilesBond[1]
			if index in smilesIndexMatches:
				continue
			
			for tempBond in templateBonds:

				tempIndex = tempBond[0]
				tempSymbol = tempBond[1]

				if symbol == tempSymbol and tempIndex not in tempIndexMatches:
					# print(smilesBond, " == ", tempBond)
					atomMatches.append([smilesBond, tempBond])
					tempIndexMatches.append(tempIndex)
					break
				# else:
					# print(smilesBond, " != ", tempBond)
			else:
				if symbol[0] not in self.BONDS:
					RgroupCount+=1

					if RgroupCount <= len(Rgroups) - 1:
						Rgroup = Rgroups[RgroupCount]
						# print(smilesBond, " == ", Rgroup)
						RgroupIndex = Rgroup[0]
						# self.removeBond(RgroupIndex,template)
						RgroupMatches.append([smilesBond, Rgroup])
						templateAtoms[RgroupIndex][0] = index
			
		if len(RgroupMatches) + len(atomMatches) != len(templateBonds):
			if template.SMILES == 'RC(=O)N(R)C(=O)R':
				print(templateBonds)
				print(len(RgroupMatches) + len(atomMatches) ," != ", len(templateBonds))
			return False

		for matchedAtoms in atomMatches:
			
			templateIndex = matchedAtoms[1][0]
			expandAtom = matchedAtoms[0]
			if template.SMILES == 'RC(=O)N(R)C(=O)R':
				print("Calling Expand Group on ", expandAtom, " from " , matchedAtoms)
			expandedGroup = self.expandGroup(expandAtom,templateIndex,template, smilesIndexMatches)
			if not expandedGroup:
				return False
		
		return True

	def getRgroups(self,atomSet):
		
		Rgroups = []

		for atom in atomSet:
			if atom[1] == 'R':
				Rgroups.append(atom)

		return Rgroups
	
	def getMainGroups(self,atomSet):
		
		mainGroups = []

		for atom in atomSet:
			if atom[1] != 'R':
				mainGroups.append(atom)

		return mainGroups
	
	def removeBond(self, index, template):
		
		for i,bonds in enumerate(template.bondData):
			for j,bond in enumerate(bonds):
				if bond[0] == index:
					del(template.bondData[i][j])
		return 0
	
	def repetitionScrub(self,functionalGroups):

		for index, group in enumerate(functionalGroups):
			groupAtoms = group.atomData
			for compareIndex,compareGroup in enumerate(functionalGroups):
				compareAtoms = compareGroup.atomData
				if all(i in compareAtoms for i in groupAtoms) and all(i in groupAtoms for i in compareAtoms) and compareIndex != index:
					del(functionalGroups[compareIndex])
		return 0
	
	def determineCyclicGroups(self,functionalGroups):
		
		for group in functionalGroups:
			groupAtoms = group.atomData
			index = 0
			while index < len(groupAtoms) - 1:
				indicies = [groupAtoms[index][0],groupAtoms[index+1][0]]
				if all(i in self.AROMATICINDICES for i in indicies):
					group.NAME = 'Aromatic' + group.NAME
					break
				elif all(i in self.CYCLICINDICES for i in indicies):
					group.NAME = 'Cyclic' + group.NAME
					break
				else:
					index+=1
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
	
	def heirarchyScrub(self,functionalGroups):
		
		index = -1
		while index != len(functionalGroups) - 1:
			
			index += 1
			group = functionalGroups[index]
			groupAtoms = group.atomData
			mainGroupAtoms = self.getMainGroups(groupAtoms)
			numMainRAtoms = len(self.getRgroups(groupAtoms))

			# print("\nFinding comparisons for ", group.NAME,"  ", groupAtoms)
			for compareIndex, compareGroup in enumerate(functionalGroups):
				
				compareAtoms = compareGroup.atomData
				mainCompareAtoms = self.getMainGroups(compareAtoms)
				numCompareRAtoms = len(self.getRgroups(compareAtoms))

				if all(i in mainCompareAtoms for i in mainGroupAtoms) and all(i in mainGroupAtoms for i in mainCompareAtoms) and compareIndex != index:
					# print("Found", mainCompareAtoms, " equal to ", mainGroupAtoms, " for ", group.NAME, " and ", compareGroup.NAME)
					if numMainRAtoms > numCompareRAtoms:
						# print("Taking ", group.NAME, " over ", compareGroup.NAME)
						del(functionalGroups[compareIndex])
						index = -1
					elif numCompareRAtoms > numMainRAtoms:
						# print("Taking ", compareGroup.NAME, " over ", group.NAME)
						del(functionalGroups[index])
						index = -1
		return 0

	