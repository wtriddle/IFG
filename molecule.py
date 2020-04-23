import re

class molecule():
	
	def __init__(self,smiles, name):
		self.SMILES = smiles
		self.NAME =  name
		self.atomRegex = re.compile(r'[a-zA-Z]')
		self.chargeRegex = re.compile(r'\+\-')
		self.chargedMol = True if len(self.chargeRegex.findall(smiles)) != 0 else False
		self.ATOMS = ['C', 'O', 'N', 'X', 'Z', 'S', 'I',
			'F', 'c', 'n', 'o', 'x', 'z', 's', 'i', 'f', 'R']
		self.BONDS = ['=','#']
		self.BRACKETS = ['[',']']
		self.CHARGES = ['+','-']
		self.NUMBERS = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
		self.LINEARSYMBOLS = self.ATOMS + self.BONDS + self.BRACKETS + self.CHARGES
		self.RINGS = self.initializeRINGS()
		self.atomData = self.initializeAtomData()
		self.bondData = self.initializeBondData()
		self.ringCount = len(self.RINGS)
		self.atomCount = len(self.atomData)

	def initializeAtomData(self):
		
		atomData = []
		atomIndex = -1
		
		for pos,symbol in enumerate(self.SMILES):
			
			if symbol in self.ATOMS:
				atomIndex+=1

				if self.SMILES[pos-1] == '[':
					chargedGroup = self.SMILES[pos-1:pos+3]
					atom = [atomIndex, chargedGroup]
				else:
					atom = [atomIndex, symbol]
				atomData.append(atom)

		return atomData
	
	def initializeBondData(self):
		
		bondData = []
		atomIndex = -1

		for pos, symbol in enumerate(self.SMILES):
			
			if symbol in self.ATOMS:
				atomIndex+=1
				bonds = []

				if self.SMILES[pos-1] == '[':
					leftBond = self.getLeftBond(LNPos = pos-1, LNIndex = atomIndex)
					rightBonds = self.getRightBonds(RNPos = pos+2, RNIndex = atomIndex)
				else:
					leftBond = self.getLeftBond(LNPos = pos, LNIndex = atomIndex)
					rightBonds = self.getRightBonds(RNPos = pos, RNIndex = atomIndex)

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
			LNIndex -=1 

			if LN == ']':
				chargedGroup = ']'

				while self.SMILES[LNPos] != '[':
					LNPos-=1
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
				scope +=1
			if scope > 0 and not parenthGroup:
				parenthGroups.append([RNIndex+1])
			if scope > 0:
				parenthGroup+=RN
			if RN == ')':
				scope-=1
			if scope == 0 and parenthGroup:
				parenthGroups[-1].append(parenthGroup)
			if RN in self.ATOMS:
				RNIndex+=1
			if scope == 0 and RN in self.NUMBERS:
				numGroup = self.numbersHandler(RNPos)
				rightBonds.append(numGroup)
			RNPos+=1
			if RNPos >= len(self.SMILES) or scope < 0:
				break
			RN = self.SMILES[RNPos]
		
		else:
			RNIndex += 1

			if RN in self.ATOMS:
				rightBonds.append([RNIndex, RN])

			elif RN in self.BONDS:
				expBond = RN
				RNPos+=1

				if self.SMILES[RNPos] == '[':
					chargedGroup = "["

					while self.SMILES[RNPos] != ']':
						RNPos+=1
						chargedGroup+=self.SMILES[RNPos]
					rightBonds.append([RNIndex,expBond + chargedGroup])

				elif self.SMILES[RNPos] in self.ATOMS:
					rightBonds.append([RNIndex, expBond + self.SMILES[RNPos]])

			elif RN == '[':
				chargedGroup = "["

				while self.SMILES[RNPos] != ']':
					RNPos+=1
					chargedGroup+=self.SMILES[RNPos]
				rightBonds.append([RNIndex,chargedGroup])
			
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
								pos+=1
								chargedGroup+=group[1][pos]
							rightBonds.append([group[0], expBond + chargedGroup])

						elif group[1][pos] in self.ATOMS:
							rightBonds.append([group[0], expBond + group[1][pos]])

					elif symbol == '[':
						chargedGroup = '['

						while group[1][pos] != ']':
							pos+=1
							chargedGroup+=group[1][pos]
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
					closePos +=1
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


	def numbersHandler(self,pos):
		
		for ring in self.RINGS:
			openInfo = ring[0]
			closeInfo = ring[1]

			if openInfo[0] == pos:
				return [closeInfo[1],closeInfo[2]]

			elif closeInfo[0] == pos:
				return [openInfo[1],openInfo[2]]
