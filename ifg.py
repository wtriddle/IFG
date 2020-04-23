from molecule import molecule
import re

class ifg(molecule):

	def __init__(self,smiles, REFCODE):
		super().__init__(smiles, REFCODE)
		self.functionalGroups = self.findFunctionalGroups()
	
	def findFunctionalGroups(self):
		functionalGroups = []
		atomIndex = -1

		for symbol in self.SMILES:
			if symbol in self.ATOMS:
				atomIndex +=1
				atom = self.atomData[atomIndex]
				groups = self.whichGroup(atom)

				for group in groups:
					functionalGroups.append(group)
				
		return functionalGroups

	def whichGroup(self, atom):
		# print("Determining groups from ", atom)
		matches = []
		symbol = atom[1]
		index = atom[0]

		for line in open('FGlist.txt','r'):
			lineInfo = re.compile(r'\S+').findall(line)
			lineInfo[0] = lineInfo[0].replace('[R]','R')
			template = molecule(lineInfo[0], lineInfo[1])
			
			if not self.chargedMol and len(self.chargeRegex.findall(template.SMILES)) != 0:
				continue

			if symbol in template.SMILES:

				expansionPoint = 0

				# R groups from temple.bondData can be removed because their expansion leads back to a main group atom from which is was expnaded from
				for tempAtom in template.atomData:
					index = tempAtom[0]
					if tempAtom[1] == 'R':
						template.bondData[index].clear()

				# Find expansion indicies or points inside template
				atomIndex = -1
				for tempSymbol in template.SMILES:
					if tempSymbol in self.ATOMS:
						atomIndex+=1
					if symbol == tempSymbol:
						expansionPoint = atomIndex
						break				

				expandGroup = self.expandGroup(atom,expansionPoint,template)
				if expandGroup:
					# if template.SMILES == 'RC(=NO)':
						# print("Matched ", template.SMILES)
					matches.append(template)
				
		return matches

	def expandGroup(self, atom, expansionPoint, template):
		
		atomMatches = []
		RgroupMatches = []
		tempIndexMatches = []

		templateBonds = template.bondData[expansionPoint]
		templateAtoms = template.atomData
		Rgroups = self.getRgroups(templateBonds)
		RgroupCount = -1
		atomIndex = atom[0]
		smilesBonds = self.bondData[atomIndex]
		if len(templateBonds) > len(smilesBonds):
			return False
		# if template.SMILES == 'RC(=NO)':
		# 	print("\nInitial Conditions:")
		# 	print("templateSmiles = ", template.SMILES)
		# 	print("smilesAtom = ", self.atomData[atomIndex])
		# 	print("smilesBonds = ",  smilesBonds)
		# 	print("expansionPoint = ", expansionPoint)
		# 	print("templateBonds = ", templateBonds)
		# 	print("templateAtoms = ", templateAtoms)
		# 	print("Rgroups = ", Rgroups)

		templateAtoms[expansionPoint][0] = atom[0]
		self.removeBond(expansionPoint,template)

		for smilesBond in smilesBonds:

			index = smilesBond[0]
			symbol = smilesBond[1]
			
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
			# print(templateBonds)
			# print(len(RgroupMatches) + len(atomMatches) ," != ", len(templateBonds))
			return False

		for matchedAtoms in atomMatches:
			
			templateIndex = matchedAtoms[1][0]
			expandAtom = matchedAtoms[0]
			self.removeBond(templateIndex,template)

			expandedGroup = self.expandGroup(expandAtom,templateIndex,template)
			if not expandedGroup:
				return False
		
		return True

	def getRgroups(self,bonds):
		
		Rgroups = []

		for bond in bonds:
			if bond[1] == 'R':
				Rgroups.append(bond)

		return Rgroups

	def removeBond(self, index, template):
		
		for i,bonds in enumerate(template.bondData):
			for j,bond in enumerate(bonds):
				if bond[0] == index:
					del(template.bondData[i][j])
		return 0