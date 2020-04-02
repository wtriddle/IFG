# Written by William Riddle from 9/23/19 - 4/1/20

import re
import sys
from dataManipulators import dataManipulators


class ifg(dataManipulators):

	# Initialization of object
	def __init__(self, inputSMILES):

		# Relevant Lists for data computation
		self.__ALCOHOLICINDICES = []
		self.__FGdata = []
		self.__ATOMS = ['C', 'O', 'N', 'X', 'Z', 'S', 'I',
			'F', 'c', 'n', 'o', 'x', 'z', 's', 'i', 'f']
		self.__CHARGES = ['+', '-']
		self.__PARENTHESIS = ['(', ')']
		self.__BRACKETS = ['[', ']']
		self.__BONDS = ['=', '#']
		self.__NUMBERS = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
		self.__NUMBERSREGEX = re.compile(r'\d')
		self.__RGROUPREGEX = re.compile(r'R')
		self.__ATOMSREGEX = re.compile(r'[a-zA-Z]')
		self.AMINOACID = True if len(re.compile(r'\[[nN]H[23]?\+\]').findall(inputSMILES)) != 0 else False

		# Computation of functional group data
		self.SMILES = self.formatSMILEScode(inputSMILES)
		self.__SMILESlength = len(self.SMILES)
		self.CYCLICINDICES = self.initializeCYCLICINDICES()
		self.RINGPOSITIONS = self.initializeRINGPOSITIONS()
		self.AROMATICINDICES = self.initializeAROMATICINDICES()
		self.FGDATAFINAL = self.identifyFunctionalGroups()

	# Identification process can be split into a two phase process
	def identifyFunctionalGroups(self):
		
		# Collection Phase of script is run first to determine all possible functional groups, placed into __FGdata
		self.collectFunctionalGroups()

		# Evaluation Phase to evalute collected groups, finalize them with accuracy and precision, and gather other various data parameters
		data = self.evaluateFunctionalGroups()

		return data
	
	# Collection Phase
	def collectFunctionalGroups(self):
			
		atomIndex = -1 # Atom index counter

		# Loop through SMILES code and call symbol handlers on every symbol to determine if groups can be build around them
		for SMILEScodePos, symbol in enumerate(self.SMILES):
			
			if symbol in self.__ATOMS:
				atomIndex += 1 # Increment index each atom symbol

				# Call atomHandler on the atom, pass in positional parameters
				directGroup = self.atomHandler(symbol, SMILEScodePos, atomIndex)
				if directGroup is not False:
					for group in directGroup:
						self.__FGdata.append(group)
				del(directGroup) # Dump the data after evaluated

				if symbol == 'O':
					# Determine if an alochol is present at a given oxygen
					self.determineAlcoholGroup(SMILEScodePos, atomIndex)
				elif symbol == 'N':
					# Determine if a primaryAmine is present at a given nitrogen
					self.determinePrimaryAmine(SMILEScodePos, atomIndex)

			elif symbol in self.__CHARGES:
				
				# Call chargeHandler on the charge and pass in positional parameters
				chargeGroup = self.chargeHandler(symbol, SMILEScodePos, atomIndex)
				if chargeGroup is not False:
					for group in chargeGroup:
						self.__FGdata.append(group)
				del(chargeGroup) # Dump the data after evaluated

			elif symbol in self.__BONDS:
				
				# Call bondHandler on the bond, pass in positional parameters
				bondGroup = self.bondHandler(symbol, SMILEScodePos, atomIndex)
				if bondGroup is not False:
					for group in bondGroup:
						self.__FGdata.append(group)
				del(bondGroup) # Dump the data after evaluated

	# Evaluation Phase
	def evaluateFunctionalGroups(self):
			
		# Scrub the FGdata list for removal of incorrect groups or accurate groups
		temp1 = self.repetitionScrub()
		temp2 = self.alcoholScrub(temp1)
		temp3 = self.heirarchyScrub(temp2)
		temp4 = temp3[:]
		for group in self.__ALCOHOLICINDICES:
			temp3.append(['Alcohol', 'ROH', [group[0]], [['R', group[0]]]])
		FGdataFinaltemp = self.containmentScrub(temp3)

		# The individual functional group lists appended to __FGdata have alias's for each other
		# To create allFGs from temp4, must copy not just temp4, but individual FG lists
		allFGstemp = []
		for innerGroup in temp4:
			allFGstemp.append(innerGroup[:])

		# Add acetlas to both FG temp lists
		acetals = self.determineAcetalGroups(allFGstemp)
		for acetal in acetals:
			allFGstemp.append(acetal)
			FGdataFinaltemp.append(acetal)

		# Count the alcohols to determine total amount and add FGInfo's of each alcohol with their carbon related index
		alcoholCount = len(self.__ALCOHOLICINDICES)

		# Finalize functional group data by determining the cyclic groups present
		FGdataFinal = self.determineCyclicGroups(FGdataFinaltemp)
		allFGs = self.determineCyclicGroups(allFGstemp)

		# Create data dictionaries to count the number of functional group occurances. Add totalAlcohols count to both of them
		FGdataFinalDict = self.createFGDataDict(FGdataFinal)
		allFGsDict = self.createFGDataDict(allFGs)
		FGdataFinalDict.update({"totalAlcohols" : alcoholCount})
		allFGsDict.update({"totalAlcohols" : alcoholCount})

		# Occurance of [NH+],[NH2+], or [NH3+] gaurnetees amino acid strucutre
		if self.AMINOACID is True:
			FGdataFinalDict.update({'AminoAcid' : 1})
		else:
			FGdataFinalDict.update({'AminoAcid' : 0})

		# Determine the number of aromatic and nonAromatic rings from RINGPOSITIONS 
		ringCount = len(self.RINGPOSITIONS)
		aromaticRingCount = nonAromaticRingCount = 0
		for ring in self.RINGPOSITIONS:
			if ring[0][1].islower(): # Lowercase symbol means aromatic ring in SMILEScode
				aromaticRingCount+=1
			if ring[0][1].isupper(): # Uppercase symbol means non-aromatic ring in SMILEScode
				nonAromaticRingCount+=1
			# Make sure the total ring count is valid, then add the data to FGdataDict

		if nonAromaticRingCount + aromaticRingCount == ringCount:
			FGdataFinalDict.update({'RingCount' : ringCount})
			FGdataFinalDict.update({'aromaticRingCount' : aromaticRingCount})
			FGdataFinalDict.update({'nonAromaticRingCount' : nonAromaticRingCount})
		else:
			print("Fatal ring count error. Exiting..")
			self.printLists()
			print(nonAromaticRingCount, " + ", aromaticRingCount, " != ", ringCount )
			print(self.SMILES)
			sys.exit()
		return (FGdataFinal, FGdataFinalDict, allFGs, allFGsDict)
		
	# Alcohol Processing
	def determineAlcoholGroup(self, SMILEScodePos, atomIndex):

		# Final atom alcohol case
		# Must be an oxygen as final positional symbol and atom in SMILEScode that is not attacthed to a ring. Must be individual
		if SMILEScodePos == self.__SMILESlength - 1:
			# May be next to a carbon, at the closure of a ring, or next to an outer group
			if self.SMILES[SMILEScodePos-1].upper() == 'C' or self.SMILES[SMILEScodePos-1] in self.__NUMBERS or self.SMILES[SMILEScodePos-1] == ')':
				# Intialze positional variables for outer and non-outer cases
				# (RRR...R(RRR..)O) case
				if self.SMILES[SMILEScodePos-1] == ')':
					numRGroups = 3
					lBlockCounter = 1
					LNPos = SMILEScodePos - 2
				else:
					numRGroups = 2
					lBlockCounter = 0
					LNPos = SMILEScodePos - 1
				# LN expansion logic to determine if an alcohol group exists
				LN = self.SMILES[LNPos]
				LNIndex = atomIndex
				while LN.upper() != 'C' or lBlockCounter > 0:
					if LN in self.__ATOMS:
						LNIndex -= 1
					if LN not in self.__NUMBERS and lBlockCounter == 0:
						break
					if LN in self.__NUMBERS and lBlockCounter == 0:
						numRGroups += 1
					if LN == '(':
						lBlockCounter -= 1
					if LN == ')':
						lBlockCounter += 1
					# Increment at bottom of loop so else statement below handles ALCOHOLICINDICES
					LNPos -= 1
					LN = self.SMILES[LNPos]
				else:
					# Append LNIndex-1 to account for stopping at the C, but not evalutating it
					self.__ALCOHOLICINDICES.append([LNIndex-1])
					self.__ALCOHOLICINDICES[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
					self.__ALCOHOLICINDICES[-1].append("FinalO")

		# First atom alcohol case
		# Check for oxygen at first position of SMILEScode, and a single carbon bound to it
		elif atomIndex == 0 and self.SMILES[SMILEScodePos+1].upper() == 'C':
			self.__ALCOHOLICINDICES.append([atomIndex+1]) # Automatically append next index if satisfied
			# Check for C1( case, indicates 4 Rgroups off of C, including O
			if self.SMILES[SMILEScodePos+2] in self.__NUMBERS and self.SMILES[SMILEScodePos+3] == '(':
				numRGroups = 4
			# Check for C1R case, indicates 3 R groups off of C, including O
			elif self.SMILES[SMILEScodePos+2] in self.__NUMBERS and self.SMILES[SMILEScodePos+3] in self.__ATOMS:
				numRGroups = 3
			# Check for C(, indicates 3 R groups off of C, including O
			elif self.SMILES[SMILEScodePos+2] == '(':
				numRGroups = 3
			# All other cases are 2 R groups
			else:
				numRGroups = 2
			self.__ALCOHOLICINDICES[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
			self.__ALCOHOLICINDICES[-1].append("First O")
		# (O) single oxygen outer groups are automatically alcoholic oxygens since they are isolated
		elif self.SMILES[SMILEScodePos-1:SMILEScodePos+2] == '(O)':
			self.__ALCOHOLICINDICES.append([atomIndex-1]) # Atom index pointing to R group bound to (O), i.e. C(O), pointing to the C
			# No loop necessary, only a couple of cases. LN variables for clarity
			LNPos = SMILEScodePos - 2 # Pointing to symbol bound directly ouside of (O) group
			LN = self.SMILES[LNPos]
			# Check for C1(O) case, indicates 4 Rgroups off of C, including O
			if LN in self.__NUMBERS:
				numRGroups = 4
			# Check for C(O)(RRR....) case, indicates 4 Rgroups off of C, including O
			elif LN in self.__ATOMS and self.SMILES[SMILEScodePos+2] == '(':
				numRGroups = 4
			# Check for C(RRR....)(O) case, indicates 4 Rgroups off of C, including O
			elif LN == ')':
				numRGroups = 4
			# An outer C(O) group indicates at least 3 always by itself if no 4 cases are satisfied
			else:
				numRGroups = 3
			self.__ALCOHOLICINDICES[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
			self.__ALCOHOLICINDICES[-1].append("Beginning (")
		# Lone oxygen closing outer group case, must be a (RRR....O) case
		elif self.SMILES[SMILEScodePos+1] == ')':
			# Initalize positional varibles of outer and non-outer cases bound leftward to O
			# (RRR...R(RRR..)O) case
			if self.SMILES[SMILEScodePos-1] == ')':
				lBlockCounter = 1
				numRGroups = 3
				LNPos = SMILEScodePos - 2
			else:
				numRGroups = 2
				lBlockCounter = 0
				LNPos = SMILEScodePos - 1
			# LN expansion logic to determine if an alcohol group exists at this particular oxygen
			LN = self.SMILES[LNPos]
			LNIndex = atomIndex
			while LN.upper() != 'C' or lBlockCounter > 0:
				if LN in self.__ATOMS:
					LNIndex -= 1
				if LN not in self.__NUMBERS and lBlockCounter == 0:
					break
				if LN in self.__NUMBERS and lBlockCounter == 0:
					numRGroups += 1
				if LN == '(':
					lBlockCounter -= 1
				if LN == ')':
					lBlockCounter += 1
				LNPos -= 1
				LN = self.SMILES[LNPos]
			else:
				# Append LNIndex-1 to account for stopping at C an not evaluating it
				self.__ALCOHOLICINDICES.append([LNIndex-1])
				self.__ALCOHOLICINDICES[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
				self.__ALCOHOLICINDICES[-1].append("Closing )")
		return 0

	# Primary Amine Processing
	def determinePrimaryAmine(self, SMILEScodePos, atomIndex):
		
		# Final atom nitrogen case
		# Must be a nitrogen as final positional symbol and atom in SMILEScode that is not attacthed to a ring. Must be individual
		if SMILEScodePos == self.__SMILESlength - 1:
			# May be next to a carbon, at the closure of a ring, or next to an outer group
			if self.SMILES[SMILEScodePos-1].upper() == 'C' or self.SMILES[SMILEScodePos-1] in self.__NUMBERS or self.SMILES[SMILEScodePos-1] == ')':
				# Intialze positional variables for outer and non-outer cases
				# (RRR...R(RRR..)N) case
				if self.SMILES[SMILEScodePos-1] == ')':
					lBlockCounter = 1
					LNPos = SMILEScodePos - 2
				else:
					lBlockCounter = 0
					LNPos = SMILEScodePos - 1
				# LN expansion logic to determine if a PrimaryAmine exists
				LN = self.SMILES[LNPos]
				LNIndex = atomIndex
				while LN.upper() != 'C' or lBlockCounter > 0:
					if LN in self.__ATOMS:
						LNIndex -= 1
					if LN not in self.__NUMBERS and lBlockCounter == 0:
						break
					if LN == '(':
						lBlockCounter -= 1
					if LN == ')':
						lBlockCounter += 1
					# Increment at bottom of loop so else statement below handles stopping at said positonal variable
					LNPos -= 1
					LN = self.SMILES[LNPos]
				else:
					# Append LNIndex-1 to account for stopping at the C, but not evalutating it
					self.__FGdata.append(['PrimaryAmine', 'RN', [LNIndex-1,atomIndex], [['R', LNIndex-1], ['N', atomIndex]]])
		# First atom nitrogen case
		# Check for nitrogen at first position of SMILEScode, and a single carbon bound to it
		elif atomIndex == 0 and self.SMILES[SMILEScodePos+1].upper() == 'C':
			self.__FGdata.append(['PrimaryAmine', 'NR', [0,1], [['N', 0], ['R', 1]]])		
		# (N) single nitrogen outer groups are automatically primary amines
		elif self.SMILES[SMILEScodePos-1:SMILEScodePos+2] == '(N)':
			self.__FGdata.append(['PrimaryAmine', 'R(N)', [atomIndex,atomIndex+1], [['R', atomIndex-1], ['N', atomIndex]]])
			# No loop necessary, only a couple of cases. LN variables for clarity
		# Lone nitrogen closing outer group case, must be a (RRR....N) case
		elif self.SMILES[SMILEScodePos+1] == ')':
			# Initalize positional varibles of outer and non-outer cases bound leftward to N
			# (RRR...R(RRR..)N) case
			if self.SMILES[SMILEScodePos-1] == ')':
				lBlockCounter = 1
				LNPos = SMILEScodePos - 2
			else:
				lBlockCounter = 0
				LNPos = SMILEScodePos - 1
			# LN expansion logic to determine if an PrimaryAmine group exists at this particular nitrogen
			LN = self.SMILES[LNPos]
			LNIndex = atomIndex
			while LN.upper() != 'C' or lBlockCounter > 0:
				if LN in self.__ATOMS:
					LNIndex -= 1
				if LN not in self.__NUMBERS and lBlockCounter == 0:
					break
				if LN == '(':
					lBlockCounter -= 1
				if LN == ')':
					lBlockCounter += 1
				LNPos -= 1
				LN = self.SMILES[LNPos]
			else:
				# Append LNIndex-1 to account for stopping at C an not evaluating it
				self.__FGdata.append(['PrimaryAmine', 'RN', [LNIndex-1,atomIndex], [['R', LNIndex-1], ['N', atomIndex]]])
		return 0

	# Acetal Processing
	def determineAcetalGroups(self, dataList):

		# Initializations
		acetalGroups = [] # Holds the acetal groups if found
		etherGroups = [] # Holds ether groups which will be looped over to determine acetals
		acetalGroupindices = [] # Tracks the index lists of each acetal group
		groupCounter = -1 # Counter for loops

		# Capture all ether groups
		for group in dataList:
			if group[0] == "Ether":
				etherGroups.append(group)

		# Loop through ether groups to determine acetals
		for group in etherGroups:
			
			groupIndices = group[2]

			# Find Hemiacetals and Hemiketals via the alcoholic carbons, which are considered acetal carbons with respect to the ether
			for acetalCarbon in self.__ALCOHOLICINDICES:
				if acetalCarbon[0] in groupIndices: # Check if the acetal carbon is part of the ether. If so, it is part of an acetal group
					numRgroups = acetalCarbon[1] # Grab the numRGroups from associated list
					if numRgroups == 3: # 3 Rgropus means Hemiacetal group
						# Create Hemiacetal group from known information
						HemiacetalIndices = []
						for index in groupIndices: # All ether indices are part of Hemiacetal indices
							HemiacetalIndices.append(index)
						HemiacetalIndices.append(acetalCarbon[0]) # As well as alcoholic, or acetal, carbon index
						distinctAcetal = not all(index in acetalGroupindices for index in HemiacetalIndices) # Check if the Hemiacetal group has already been found
						if distinctAcetal is True: # If not already found, then it is a new Hemiacetal group. Add its contents to acetalGroups
							if HemiacetalIndices not in acetalGroupindices:
								acetalGroups.append(['Hemiacetal', '', HemiacetalIndices ,[[]]])
								for index in HemiacetalIndices:
									acetalGroupindices.append(index)
								acetalGroupindices.append(HemiacetalIndices)
						break
					elif numRgroups == 4: # 4 Rgroups means Hemiketal group
						# Create Hemiketal group from known information
						HemiketalIndices = []
						for index in groupIndices: # All ether indices are part of Hemiketal indices
							HemiketalIndices.append(index)
						HemiketalIndices.append(acetalCarbon[0]) # As well as alcoholic, or acetal, carbon index
						distinctAcetal = not all(index in acetalGroupindices for index in HemiketalIndices) # Check if the Hemiketal group has already been found
						if distinctAcetal is True: # If not already found, then it is a new Hemiketal group. Add its contents to acetalGroups
							if HemiketalIndices not in acetalGroupindices:
								acetalGroups.append(['Hemiketal', '', HemiketalIndices ,[[]]])
								for index in HemiketalIndices:
									acetalGroupindices.append(index)
								acetalGroupindices.append(HemiketalIndices)
						break
					else: # 2 Rgroups means standalone alcohol. Go to next ether group by breaking the ALCOHOLICINDICES loop
						break

			# Acetal Loop
			# Acetals require two ethers with proper crossover, so a loop over etherGroups is used
			for compareGroupCounter, compareGroup in enumerate(etherGroups):

				compareIndices = compareGroup[2]

				# If the ethers are the same, skip to next one
				if compareGroupCounter == groupCounter:
					continue

				# If the oxygens of an ether are the same as another, then so are the R groups, and thus they are the same ether group. 
				# Check for identical oxygens to check for identical ether groups
				identicalOxygens = all(index in compareIndices for index in groupIndices)
				if identicalOxygens is True:
					continue

				# If the oxygens are distinct, then the groups are not identical, but they either have 1 or 0 similar R groups possible
				# Thus, if there is any crossover without identical oxygens, then there is an acetal group present
				acetalGroup = any(index in compareIndices for index in groupIndices)
				if acetalGroup is True:
					# Combine ether indices
					acetalIndices = []
					for index in groupIndices:
						acetalIndices.append(index)
					for index in compareIndices:
						if index not in acetalIndices:
							acetalIndices.append(index)
					# Dictionary and template are omitted. Acetals do not have distinct templates in the SMILEScodes,
					# And the templates are not necessary to determine cyclic or non-cyclic distinctions, only indices
					# Add a distinctly found acetal group to indices
					distinctAcetal = not all(index in acetalGroupindices for index in acetalIndices) # Check if the acetal group has already been found
					if distinctAcetal is True:
						acetal = ["Acetal", "", acetalIndices, []]
						for index in acetalIndices:
							acetalGroupindices.append(index)
						acetalGroups.append(acetal)
					break
		return acetalGroups # Returns list of acetal groups

	# Create dynamic strucutre from SMILES code
	def initializeRINGPOSITIONS(self):

		# Initializations
		atomIndex = -1 # Index counter
		evaluatedNumbers = [] # List that tracks the string positions of the numbers whose ring junction information has already been determined
		RINGPOSITIONS = []

		# Main Loop
		for SMILEScodePos, symbol in enumerate(self.SMILES):

			# Set two properties which distinguish if a number is an outer group, i.e. (R), or a linear group, i.e. R.
			# See documentation for explanation of this SMILEScode analysis
			openingProperty = closingProperty = ""
			
			# Keep track of atom positional data for when a number is found
			if symbol in self.__ATOMS:
				atomIndex += 1
				# Atom found before a new, unevalutaed number is the "opening" atom of that ring.
				openingAtom = symbol
				openingAtomPos = SMILEScodePos

			# When an unevalutaed number, or essentially a new ring, is encountered, evalute it
			if symbol in self.__NUMBERS and SMILEScodePos not in evaluatedNumbers:
				chargeGroup = "" # Holds charged atoms bracket group, i.e. [N+] if a chargeGroup opens or closes a ring. Reset to an empty string for subsequent charge groups

				# The openingAtomPos variable points to the charged atom inside the brackets for the following scenario.
				# A closing bracket must precede the position of the number instead of an atom if a charge group opens the ring.
				if self.SMILES[SMILEScodePos-1] == ']' and SMILEScodePos >= openingAtomPos + 2:
					specialCounter = SMILEScodePos - 1 # Initialized at position of left bracket
					specialSymbol = self.SMILES[specialCounter]
					while specialSymbol != '[':
						chargeGroup += specialSymbol
						specialCounter -= 1
						# If the counter goes out of the scope of the SMILEScode, exit the script
						if specialCounter == -1:
							print("Fatal error in the creation of RINGPOSITIONS opening charge group on the smilescode ", self.SMILES)
							sys.exit()
						specialSymbol = self.SMILES[specialCounter]
					chargeGroup += '[' # Loop ceases on '[', so add it to the chargeGroup
					chargeGroup = chargeGroup[::-1] # Reverse the string because it was appended leftward
					# Reassign openingAtom and openingAtomPos to the charge group and special counter
					# so that the data can be appended to RINGPOSITIONS with a single statement instead of a lengthy if-else chargeGroup statement
					openingAtom = chargeGroup
					openingAtomPos = specialCounter

				# Initialize subLoop variables
				correlatingPos = SMILEScodePos
				correlatingIndex = atomIndex

				# SubLoop to find correlating atom at ring junction. Begin one symbol after the number
				for subSymbol in self.SMILES[SMILEScodePos+1:self.__SMILESlength]:
					correlatingPos += 1 # Track string position

					# Keep track of atom positional data for when the correlating atom is found
					if subSymbol in self.__ATOMS:
						correlatingAtom = subSymbol
						correlatingIndex += 1
						correlatingAtomPos = correlatingPos

					# If the same number which began the sub loop is found, then the ring closure has been found
					# The correlatingAtom positional data and opening atom data will be paired together in RINGPOSITIONS
					if subSymbol == symbol:
						# If the ring opens with a charge, it should close with a charge. Re-assign positional variables
						# This special case is evaluated in the same manner as in the opening case
						if chargeGroup != "":
							correlatingChargeGroup = ""
							specialCounter = correlatingPos -1
							specialSymbol = self.SMILES[specialCounter]
							while specialSymbol != '[':
								correlatingChargeGroup += specialSymbol
								specialCounter -= 1
								# If the counter goes out of the scope of the SMILEScode, exit the script
								if specialCounter == -1:
									print("Fatal error in the creation of RINGPOSITIONS correlating charge group on the smilescode ", self.SMILES)
									sys.exit()
								specialSymbol = self.SMILES[specialCounter]
							correlatingChargeGroup = correlatingChargeGroup[::-1] # Flip because of leftward expansion

							# Charge group must be a length of at most 4. Anything larger suggests an error, so take the other atom first
							if len(correlatingChargeGroup) <= 4:
								# Reset the variables in this special case
								correlatingAtom = correlatingChargeGroup
								correlatingAtomPos = specialCounter

						# Determine linear or outer distinctions for each atom
						# See documentation for greater detail of group distinction

						# Opening property distinction
						RNPos = openingAtomPos
						RN = self.SMILES[RNPos]
						while openingProperty == "":
							RNPos+=1
							if RNPos >= self.__SMILESlength - 1:
								openingProperty = "Linear"
								break
							RN = self.SMILES[RNPos]
							if RN in self.__ATOMS or RN in self.__BONDS:
								openingProperty = "Outer"
								break
							if RN == ')':
								openingProperty = "Linear"
								break

						# Correlating or closing property distinction
						RNPos = correlatingPos
						RN = self.SMILES[RNPos]
						while closingProperty == "":
							RNPos+=1
							if RNPos >= self.__SMILESlength - 1:
								closingProperty = "Linear"
								break
							RN = self.SMILES[RNPos]
							if RN in self.__ATOMS or RN in self.__BONDS:
								closingProperty = "Outer"
								break
							if RN == ')':
								closingProperty = "Linear"
								break

						# Append the opening and closing (or correlating) information into a single list.
						# Each list within RINGPOSITIONS indicates one set of ring information
						openingInfo = [SMILEScodePos, openingAtom, atomIndex, openingAtomPos, openingProperty]
						closingInfo = [correlatingPos, correlatingAtom, correlatingIndex, correlatingAtomPos, closingProperty]
						RINGPOSITIONS.append([openingInfo, closingInfo])

						# Add numbers to evaluted list
						evaluatedNumbers.append(SMILEScodePos)
						evaluatedNumbers.append(correlatingPos)
						break

					# Error caught if correlatingPos runs out of the scope of the SMILEScode
					if correlatingPos == self.__SMILESlength - 1:
						print("There has been an error in the creation of the RINGPOSITIONS global \n")
						print("The SMILEScode code ", self.SMILES, " had an error")
						sys.exit()
		return RINGPOSITIONS

	# DLA to SAR conversion + hydrogen removal for initializing the self.SMILES attribute
	def formatSMILEScode(self, SMILEScode):
		SMILEScodePos = -1
		reFormatted = ""
		for symbol in SMILEScode:
			SMILEScodePos += 1
			if symbol == '[':
				startBracketPos = SMILEScodePos
			if symbol == 'H':
				cutPos = SMILEScodePos
				while SMILEScode[cutPos] != ']':
					cutPos += 1
				reFormatted = SMILEScode[0:startBracketPos] + SMILEScode[startBracketPos+1] + SMILEScode[cutPos+1:len(SMILEScode)]
				reFormatted = self.formatSMILEScode(reFormatted)
				break
		if reFormatted == "": # If no formatted occured, DLA to SAR the original SMILEScode input
			return self.DLAtoSARconversion(SMILEScode)
		else: # Otherwise, give reFormatted to DLA to SAR converter
			return self.DLAtoSARconversion(reFormatted)

		# DLA to SAR conversion

	# Double Lettered Atom to Single Atom Representation converter
	def DLAtoSARconversion(self, template):
		# conversionLegend
		# BR --> X
		# CL --> Z
		templatePos = -1
		reformattedTemplate = ""
		while templatePos != len(template) - 1:
			templatePos += 1
			# If symbol is doubled lettered, change it to a single letter representation
			if template[templatePos:templatePos+2] == "Br":
				reformattedTemplate += "X"
				templatePos += 1
			elif template[templatePos:templatePos+2] == "Cl":
				reformattedTemplate += "Z"
				templatePos += 1
			else:
				reformattedTemplate += template[templatePos]
		return reformattedTemplate

	# Create the CYCLICINDICES list
	def initializeCYCLICINDICES(self):

		# Tracks the position of the numbers inside of the SMILEScode to not repeat ring analysis when closing ring is found
		evaluatedNumbers = []
		CYCLICINDICES = []

		atomIndex = -1 # Atom counter
		for SMILEScodePos, symbol in enumerate(self.SMILES):

			if symbol in self.__ATOMS:
				atomIndex += 1

			# If a new ring has been found in the SMILES code via a number 
			if symbol in self.__NUMBERS and SMILEScodePos not in evaluatedNumbers:
				# print("Found ", symbol, " at position, ", SMILEScodePos)
				# print("\n")

				outerIndices = []
				rBlockCounter = 0 # Counter to ensure same scope of indices is appended to CYCLICINDICES
				evaluatedNumbers.append(SMILEScodePos) # Index directly next to number is always in ring
				if atomIndex not in CYCLICINDICES:
					CYCLICINDICES.append(atomIndex)
				RNPos = SMILEScodePos + 1
				RNindex = atomIndex
				# print(SMILEScode)
				RN = self.SMILES[RNPos]

				# Symbol is equal to number which opened the ring. Loop untl the same number is encountered, held in the RN positional variable
				while RN != symbol:
					# print("RN = ", RN)
					# print("rBlockCounter = ", rBlockCounter)
					if RN == '(':
						rBlockCounter += 1
					if RN == ')':
						outerIndices[rBlockCounter-1].clear()
						rBlockCounter -= 1
					if rBlockCounter == 0:
						outerIndices.clear() # Refresh every scope reset
					if RN in self.__ATOMS:
						RNindex+=1
						if rBlockCounter == 0 and RNindex not in CYCLICINDICES:
							CYCLICINDICES.append(RNindex)
						elif rBlockCounter != 0: # If an atom is within a parenthesis scope
							if len(outerIndices) < rBlockCounter: # If a new nested parenthesis or deeper parenthesis is found, i.e. in between a parenthesis group (...()..,)
								outerIndices.append([RNindex]) # Create a new list tracking the inner indices of that nested parenthesis group
							else: # If the parenthesis scope is not deeper, then access the current scope and append the atom index to it
								outerIndices[rBlockCounter-1].append(RNindex)
					RNPos+=1
					RN = self.SMILES[RNPos]
				else:
					# print("Correlating number", RN, " found at ", RNPos)
					evaluatedNumbers.append(RNPos) # Index directly next to number is always part of ring
					while rBlockCounter != 0: # If ceased within a parenthesis scope
						rBlockCounter -= 1
						for index in outerIndices[rBlockCounter]: # Add all indices within that scope to CYCLICindices
							if index not in CYCLICINDICES:
								CYCLICINDICES.append(index)
		return CYCLICINDICES
		
	# Create the AROMATICINDICES list
	def initializeAROMATICINDICES(self):
		AROMATICINDICES = []
		# Distinguish between aromatic and non aromatic indices
		atomIndex = -1
		for symbol in self.SMILES:
			if symbol in self.__ATOMS:
				atomIndex+=1
				# Lower case atoms are always aromatic
				if symbol.islower() and atomIndex not in AROMATICINDICES:
					AROMATICINDICES.append(atomIndex)
		return AROMATICINDICES

	# Scrubs self.__FGdata for repeated groups
	def repetitionScrub(self):

		# Initializations
		formattedData = [] # Captures the correct groups from after list is evaluted

		# Removes repetitions of certain FG's in dataList
		# The result is placed into formattedData
		for group in self.__FGdata:
			inFormattedData = self.isGroupInList(formattedData, group) # Check for group equivalnce in dataList
			if inFormattedData is True:
				continue
			else:
				formattedData.append(group)

		return formattedData

	# Scrubs the FG dataList for any alcohol containing groups which in fact do not point to, or contain an, alcohol. i.e Carboxlyic acid with no alcohol
	def alcoholScrub(self, dataList):

		formattedData = dataList[:] # Holds the corrected list after alcohol evalution
		groupCounter = - 1 # Counter

		while groupCounter <= len(formattedData):
			groupCounter += 1 # Track index of group

			# If the loop limit has been reached, then exit the loop
			if groupCounter == len(formattedData):
				break

			# Grab the group that groupCounter points to
			group = formattedData[groupCounter]
			groupIndices = group[2]

			# If the group contians an alcohol
			if len(re.compile(r'acid|oxide|oxime').findall(group[0])) != 0:
				# Loop through each indivudal alcohol in ALCOHOLICINDICIES to find if an individual alcohol group in the SMILEScode correlates to an alcohol in the functional group
				# i.e. a Carboxylic acid is identified, prove that there is an alcohol in the group
				# This is necessary because alcohol groups are implied and not distinguished between inside handler functions
				validGroup = False # False until proven True
				for alcohol in self.__ALCOHOLICINDICES:
					validGroup = any(index in groupIndices for index in alcohol) # If the alcohol containing group has an ALCOHOLICINDEX, then it is valid
					if validGroup is True: # If the group is valid, break the loop and continue evalution for other groups
						break
				# If the group which required an alocholic index did not have an alocholic index, then the group is False
				# Discard its contents from formattedData and reset the loop
				if validGroup is False:
					# print("REMOVING NONALCOHOL GROUP", formattedData[groupCounter])
					formattedData.pop(groupCounter)
					groupCounter = -1
					continue
		return formattedData

	# Scrubs FG dataList for heirarchichal functional groups that may be removed, i.e. T. Amine -> S. Amine -> P. Amine hierarhcy analysis
	def heirarchyScrub(self, dataList):

		formattedData = dataList[:] # Holds data after hierarcical scrub
		groupCounter = -1 # Counter

		while groupCounter <= len(formattedData):
			groupCounter += 1 # Track index of group

			# If the loop limit has been reached, then exit the loop
			if groupCounter == len(formattedData):
				break
			# Grab the group that groupCounter points to
			group = formattedData[groupCounter]
			groupIndices = group[2]

			# Loop over the hierarchy lines in FGheirarchy file that shows the proper heirarchy of functional groups
			for line in open('FGheirarchy.txt', 'r'):

				# Collect heirarchy information
				lineInfo = re.compile(r'\S+').findall(line)
				heirarchy = lineInfo[0].split(":")

				# Check if the current group is in that heirarchy
				if group[0] in heirarchy:

					# If a group was found to be part of a hierarchy, try to find another group in the same heirarchy
					for compareGroupCounter, compareGroup in enumerate(formattedData):

						# Find two distinct groups in the same heirarchy
						if compareGroup[0] in heirarchy and compareGroupCounter != groupCounter:
							compareIndices = compareGroup[2]

							# If group is a higher order functional group
							if heirarchy.index(group[0]) >= heirarchy.index(compareGroup[0]):
								# Check if the groups overlap inside of the SMILES
								fullContainment = all(index in groupIndices for index in compareIndices)
								if fullContainment is True:
									# print("REMOVING LOWER HIERARCHY", formattedData[compareGroupCounter])
									formattedData.pop(compareGroupCounter)
									groupCounter = -1
									break

							# Or if compareGroup is a higher order functional group
							elif heirarchy.index(group[0]) <= heirarchy.index(compareGroup[0]):
									# Check if the groups overlap inside of the SMILES
								fullContainment = all(index in compareIndices for index in groupIndices)
								if fullContainment is True:
									# print("REMOVING LOWER HIERARCHY", formattedData[groupCounter])
									formattedData.pop(groupCounter)
									groupCounter = -1
									break

				else: # If the group is not in the hierarchy, go to the next line
					continue
		return formattedData

	# Scrubs FG dataList for full containment of groups. i.e. Ester takes precedence over ketones and ethers
	def containmentScrub(self, dataList):

		# Initializations
		formattedData = dataList[:] # Copy input dataList
		groupCounter = -1 # Counter

		# Continue until formattedData is looped through without any revisions
		# Removal of any groups causes a reset of the loop to the beginning for another comparison loop
		groupCounter = -1
		while groupCounter <= len(formattedData):

			groupCounter += 1 # Track index of group

			# If the loop limit has been reached, then exit the loop
			if groupCounter == len(formattedData):
				break

			# Grab the group that groupCounter points to
			group = formattedData[groupCounter]
			groupTemplate = group[1]
			groupIndices = group[2]
			groupDict = group[3]
			compareGroupCounter = -1 # Set compareGroupCounter for each new group

			# Determine number of Maingroup atoms in group. Real letters - R's
			numMainGroupAtoms = len(groupIndices) - len(self.__RGROUPREGEX.findall(groupTemplate))

			# print("Determing Containment for ", group)
			# Compare this group with all other groups found in SMILES for containment
			for compareGroup in formattedData:

				compareGroupCounter += 1 # Track index
				compareTemplate = compareGroup[1]
				compareIndices = compareGroup[2]
				compareDict = compareGroup[3]

				# Determine number of Maingroup atoms in compareGroup. Real letters - R's
				numMainCompareAtoms = len(compareIndices) - len(self.__RGROUPREGEX.findall(compareTemplate))

				# Determine if all of the group atoms are contained within the compareGroup
				# i.e. check for full group containment, like a Ketone contained within an Ester
				# All of its indices must be within a larger groups indices
				fullGroupContainment = all(index in compareIndices for index in groupIndices)

				# If group is fully contained within compareGroup, that is every atom is in comparegroup, remove the group and restart loop without that group
				if fullGroupContainment is True and len(groupIndices) < len(compareIndices):
					formattedData.pop(groupCounter)
					groupCounter = -1
					break

				# Similalry, check if comapreGroup satisfies this condition
				fullCompareGroupContainment = all(index in groupIndices for index in compareIndices)

				# If compareGroup is fully contained within a group, that is every atom is in group, remove the group and restart loop without that group
				if fullCompareGroupContainment is True and len(compareIndices) < len(groupIndices):
					formattedData.pop(compareGroupCounter)
					groupCounter = -1
					break

				# Containment also occurs if all maingroup atoms of one group are apart of the other
				# Determine mainGroupCrossover if there is no full containments
				mainGroupCrossover = 0
				for atom in groupDict:
					for compareAtom in compareDict:
						# If the atoms are symbol equivalent, and both are main group atoms
						if atom[0] == compareAtom[0] and atom[0] != 'R' and compareAtom[0] != 'R':
							# And if the atoms are the same index within the SMILEScode
							if int(atom[1]) == int(compareAtom[1]):
								# Then atoms must be of the exact same maingroup in the structue
								mainGroupCrossover += 1
				# If there is main group no relation, then group is independent of compareGroup.
				if mainGroupCrossover == 0:
					continue

				# If there is full mainGroupCrossover between the two groups, i.e. they have identical mainGroupAtoms
				if mainGroupCrossover == numMainCompareAtoms == numMainGroupAtoms:

					# If compareIndices has more R groups, take it over group
					if len(groupIndices) < len(compareIndices):
						formattedData.pop(groupCounter)
						groupCounter = -1
						break

					# If groupIndices has more R groups, take it over compareGroup
					if len(groupIndices) > len(compareIndices):
						formattedData.pop(compareGroupCounter)
						groupCounter = -1
						break

				# If all main group atoms in compareGroup are contained within group and both groups are the same length, but group has more main group atoms, then remove compareGroup
				if mainGroupCrossover == numMainCompareAtoms and numMainGroupAtoms > numMainCompareAtoms:
					formattedData.pop(compareGroupCounter)
					groupCounter = -1
					break

				# If all main group atoms in group are contained within compareGroup and both groups are the same length, but compareGroup has more main group atoms, then remove group
				if mainGroupCrossover == numMainGroupAtoms and numMainCompareAtoms > numMainGroupAtoms:
					formattedData.pop(groupCounter)
					groupCounter = -1
					break

		# Returned as FGdataFinal in IFG head function
		return formattedData # Returns formattedData

	# Determines if two functional groups are exactly identical
	def isGroupInList(self, dataList, group):

		# If there is a group in dataList, comparisons can be made
		if dataList:
			for compareGroup in dataList:

				# Comparative Variables
				compareGroupTemplate = compareGroup[1]
				compareGroupindices = compareGroup[2]
				compareGroupDict = compareGroup[3]
				numMainCompareAtoms = len(compareGroupindices) - len(self.__RGROUPREGEX.findall(compareGroupTemplate))

				# Group Variables
				groupTemplate = group[1]
				groupindices = group[2]
				groupDict = group[3]
				numMainGroupAtoms = len(groupindices) - len(self.__RGROUPREGEX.findall(groupTemplate))

				# Check for equivalent number of atoms
				if len(compareGroupindices) != len(groupindices):
					continue
				else:
					mainGroupCrossover = 0
					for atom in groupDict:
						for compareAtom in compareGroupDict:
							# Check if atoms are maingroup equivalent, and if they are in the same positional index.
							# If so, the atoms in both functionalGroup scopes are "identical" add them to mainGroupCrossover
							if atom[0] == compareAtom[0] and atom[0] != 'R':
								if int(atom[1]) == int(compareAtom[1]):
									mainGroupCrossover += 1

					# If there is full maingroups crossover and the templates are the same length
					if mainGroupCrossover == numMainCompareAtoms == numMainGroupAtoms:
							# Then the exact group is already contained within dataList
							return True
					else: # If not equivalent lengths, continue
						continue
		return False # Returns True if group is in the list, otherwise False

	# Determines the Cylic/Aromatic distinctions of functional groups within a SMILES structure
	def determineCyclicGroups(self, dataList):
		
		formattedData = []
		copiedList = dataList[:] # Copy list so original input list is not tampered with
		for group in copiedList:
			
			groupIndices = group[2]

			# If all indices of a group are Cyclic/Aromatic, then the group is automatically Cyclic/Aromatic
			isAromatic = all(index in self.AROMATICINDICES for index in groupIndices)
			isCyclic = all(index in self.CYCLICINDICES for index in groupIndices)

			# Check alcohol groups first, they are special case
			if group[0] == 'Alcohol':
				# Check for Aromatic first, since aromatic group implies cyclic group
				if groupIndices[0] in self.AROMATICINDICES:
					group[0] = "Aromatic" + group[0]
				elif groupIndices[0] in self.CYCLICINDICES:
					group[0] = "Cyclic" + group[0]

			# Check for Aromatic first, since aromatic group implies cyclic group
			elif isAromatic is True:
				group[0] = "Aromatic" + group[0]

			# Check for case of all indices being cyclic, automatically a cyclic group
			elif isCyclic is True:
				group[0] = "Cyclic" + group[0]

			# Check for case when at least 2 consecutive indices are cyclic, then group is cyclic
			else:
				sampleIndices = []
				indexCounter = -1
				while isCyclic is False:
					indexCounter += 1
					# Checking one ahead, and 0 index of lists. So break at length - 2
					if indexCounter == len(groupIndices) - 2:
						break
					# Add 2 consecutive indices to sampleIndices to see if they are both cyclic. 
					sampleIndices.append(groupIndices[indexCounter])
					sampleIndices.append(groupIndices[indexCounter+1])
					isCyclic = all(index in self.CYCLICINDICES for index in sampleIndices) # Check if the 2 consecutive indices were cyclic
					sampleIndices.clear() # Refresh for a new set each increment of indexCounter
					# If two consecutive indices were cyclic, then the group is cyclic
					if isCyclic is True:
						group[0] = "Cyclic" + group[0]
						break
			# Append the newly distinguished group to formattedData, and all other unchanged groups as well. 
			# Resultant list is same functional group list, but with Cylic/Aromatic distinctions were necessary
			formattedData.append(group)
		return formattedData

	# Print data Lists associated with this particular SMILES code
	def printLists(self):
		print("FGdata = ", self.__FGdata)
		print("ALCOHOLICINDICES = ", self.__ALCOHOLICINDICES)
		print("CYCLICINDICES = ", self.CYCLICINDICES)
		print("AROMATICINDICES = ", self.AROMATICINDICES)
		print("RINGPOSITIONS = ", self.RINGPOSITIONS)
	
	# Connects two atoms by a double or triple bonds only. Called when = or # is encountered
	def bondHandler(self, bondSymbol, bondPostion, atomIndex):

		# Initialize positional variables
		LNPos = RNPos = bondPostion # Initialized at the poistion of the bond
		LNindex = RNindex = atomIndex # Initialized at the position of the previous atom to the left of where the bond appears
		LN = RN = bondSymbol # Initialized as bondSymbol

		# Special case which occurs commonly as R(=R)
		if self.SMILES[LNPos-1] == '(' and self.SMILES[LNPos-2] in self.__ATOMS and self.SMILES[RNPos+1] in self.__ATOMS and self.SMILES[RNPos+2] == ')':
			atomindices = [atomIndex, atomIndex+1]
			bondGroup = self.SMILES[LNPos-2] + "(=" + self.SMILES[RNPos+1] + ')'
			FGinfo = self.whichGroup(LNPos-2, RNPos+2, bondGroup, atomindices)
			# Bond has been handled, so return once evaluted
			if FGinfo is not False:
				return FGinfo
			else:
				return False


		# LN loop

		# Initializations
		lBlockCounter = 0 # Tracks number of blocking groups
		lOuterParenthGroup = "" # Tracks parenthesis group blocking LN atom in bond i.e C(R)=O
		lOuterNumGroup = "" # Tracks a number group blocking LN atom in bond i.e. C1=O
		lOuterindices = [] # Trakcs the indicies inside parenthesis blocking group
		LNindex+=1 # Increment index to allow loop to account for offset of intialized atomIndex, which is already the lefthand index.
		# Therefore, must increment so that when an atom is encounter in the loop it actually points to the correct atom

		# Find the first non blocking group atom connected to the bond
		while LN not in self.__ATOMS or lBlockCounter > 0:
			LNPos -= 1
			LN = self.SMILES[LNPos]
			if LN in self.__CHARGES and lBlockCounter == 0: # Let charge handler deal with bracket charge groups
				return False
			if LN == ')':
				lBlockCounter += 1
			if LN == '(' and lBlockCounter != 0:
				lBlockCounter -= 1
			if LN in self.__NUMBERS and lBlockCounter == 0:
				lOuterNumGroup = self.numbersHandler(LNPos)
			if lBlockCounter != 0:
				lOuterParenthGroup += LN
			if LN in self.__ATOMS:
				LNindex -= 1
			if lBlockCounter != 0 and LN in self.__ATOMS:
				lOuterindices.append(LNindex)

		# RN loop
		# The atom next to the bond symbol is always the atom involved in the bond, no blocking groups ever occur
		while RN not in self.__ATOMS:
			RNPos += 1
			RN = self.SMILES[RNPos]
			if RN in self.__BRACKETS: # Let charge handler deal with bracket charge groups
				return False
			if RN in self.__ATOMS:
				RNindex += 1

		# lOuterParenthGroup evaluation to see if it should be included in bondGroup
		lOuterParenthGroup = lOuterParenthGroup[::-1] # Flip because appended from the left, so group appears backwards
		atomCount = len(self.__ATOMSREGEX.findall(lOuterParenthGroup))
		numCount = len(self.__NUMBERSREGEX.findall(lOuterParenthGroup))
		# If the outer group is cumbersome, don't include its large contents
		# Only looking for simple cases like C(O)=O or C(N1)=O, not C(R2R3R.....)=O
		if atomCount > 1 or numCount != 0:
			lOuterParenthGroup = ""

		# Bond Group creation according to the specific case satisfied
		if lOuterNumGroup and lOuterParenthGroup:
			atomindices = [LNindex, lOuterNumGroup[1], lOuterindices[0], RNindex]
			bondGroup = LN  + '(' +lOuterNumGroup[0] + ')' + lOuterParenthGroup + bondSymbol + RN
		elif lOuterNumGroup and not lOuterParenthGroup:
			atomindices = [LNindex, lOuterNumGroup[1], RNindex]
			bondGroup = LN + '(' + lOuterNumGroup[0] + ')' + bondSymbol + RN
		elif lOuterParenthGroup and not lOuterNumGroup:
			atomindices = [LNindex, lOuterindices[0], RNindex]
			bondGroup = LN + lOuterParenthGroup + bondSymbol + RN
		else:
			atomindices = [LNindex, RNindex]
			bondGroup = LN + bondSymbol + RN

		# Determine the group and return the info
		FGinfo = self.whichGroup(LNPos, RNPos, bondGroup, atomindices)
		return FGinfo # Returns False for no matches, otherwise the whichGroup info is returned based on bond expansion

	# Connects two atoms by a single bond, but may also handle double or triple when necessary. Called when any letter is encountered
	def atomHandler(self, atomSymbol, atomPosition, atomIndex):

		# Initialize Variables, next to positional call location from ifg
		LNPos = atomPosition - 1
		RNPos = atomPosition + 1
		RNindex = LNindex = atomIndex
		RNbond = LNbond = ""
		FGinfo = []

		# RN LOOP

		# Validate scope before expansion
		if RNPos < self.__SMILESlength - 1:
			# Initializations
			RN = self.SMILES[RNPos]
			# If atom handler called from within a charged and bracketed symbol, return False
			if RN in self.__CHARGES: # Let charge handler deal with charged groups
				return False
			rBlockCounter = 0 # Tracks the current number of outer parenthesis groups blocking the neighbor
			rOuterGroups = 0 # Tracks total amount of indivudal parenthesis groups blocking the neighbor
			rInnerParenthGroup = "("
			rInnerParenthindices = [] # Tracks parenthesis group indices
			rInnerParenthGroupPositions = [] # Tracks parenthesis group position in SMILEScode
			numGroup = "" # Holds number group if found

			# Loop continues until an atom within the same grouping layer as atomSymbol is found to the right
			while rBlockCounter > 0 or RN not in self.__ATOMS:
				# If RNPos reaches final position in loop without stopping, no RN atom exist. Cease the loop
				if RNPos == self.__SMILESlength - 1:
					RN = ""
					break
				if RN in self.__ATOMS:
					RNindex += 1 # Increment atom index even inside parenthesis groups
				# If scope of RN loop is within a parenthesis, capture the information inside of it with respect to atomSymbol
				if rBlockCounter != 0:
					rInnerParenthGroup += RN
					rInnerParenthGroupPositions.append(RNPos)
					if RN in self.__ATOMS:
						rInnerParenthindices.append(RNindex)
					if RN in self.__NUMBERS:
						numGroup = self.numbersHandler(RNPos)
						rInnerParenthindices.append(numGroup[1])
						if numGroup[3] == "Linear":
							rInnerParenthGroup += numGroup[0]
						elif numGroup[3] == "Outer":
							outerGroup = '(' + numGroup[0] + ')'
							rInnerParenthGroup += outerGroup
				# Only add bond within the scope of where atom was initially called from
				if RN in self.__BONDS and rBlockCounter == 0:
					RNbond = RN
				if RN == '(':
					# If a new group has been found starting at the scope of the initial atom, increment rOuterGroups count
					if rBlockCounter == 0:
						rOuterGroups += 1
					rBlockCounter += 1
				if RN == ')':
					rBlockCounter -= 1
				# If RN on the edge of a parenthesis, no RN atom exists, cease the loop
				if rBlockCounter < 0:
					RN = ""
					break
				RNPos += 1
				RN = self.SMILES[RNPos]
			else:
				RNindex += 1 # Increment index once atom found because not evaluted, but stopped on instead

			atomCount = len(self.__ATOMSREGEX.findall(rInnerParenthGroup)) # Find number of atoms inside of outer parenthesis group

			# R(RRR...)R case, one outer group and many inside atoms
			# Creates two expansion possibilities, see documentation on RN for explanation
			if atomCount >= 2 and rOuterGroups == 1:
				temp = ""
				for symbol in rInnerParenthGroup:
					# Strip parenthesis in this >2 atom group, treat as linear group
					if symbol != '(':
						temp += symbol

					# Break after first atom has been added
					if symbol in self.__ATOMS:
						break

				# Create Inner and Outer expansion group possibilities
				outerExpandGroup = atomSymbol + '(' + temp + ')' # first atom (temp) is inner, expand from outside larger parenthesis group
				innerExpandGroup = atomSymbol + '(' + RN + ')' + temp # RN is outer, expand from inside larger parenthesis group

				# Expands outside parenthesis group, with temp as parenthesied outer group attacthed to group argument
				OuterExpansion = self.whichGroup(atomPosition, rInnerParenthGroupPositions[-1], outerExpandGroup, [atomIndex, rInnerParenthindices[0]])
				if OuterExpansion is not False:
					for group in OuterExpansion:
						FGinfo.append(group)
				del(OuterExpansion)

				# Expands inside parenthesis group, with RN as parenthesied outer group attacthed to group argument
				InnerExpansion = self.whichGroup(atomPosition, rInnerParenthGroupPositions[0], innerExpandGroup, [atomIndex, RNindex, rInnerParenthindices[0]])
				if InnerExpansion is not False:
					for group in InnerExpansion:
						FGinfo.append(group)
				del(InnerExpansion)
			elif rOuterGroups > 2:
				self.printLists()
				print("Fatal error, more than two outer groups in atomHandler during RN expansion...")
				print(self.SMILES)
				sys.exit()

			# Combine RN info and determine whichGroup. If RN is empty, do not run this portion of code
			if RN != "":
				RNgroup = atomSymbol + RNbond + RN
				RNindices = [atomIndex, RNindex]
				info = self.whichGroup(atomPosition, RNPos, RNgroup, RNindices)
				if info is not False:
					for group in info:
						FGinfo.append(group)

		# LN LOOP

		# Validate SMILEScode scope
		if LNPos > 0:
			# Initializations
			LN = self.SMILES[LNPos]
			lOuterNumGroup = "" # Holds number group if encountered
			lBlockCounter = 0 # Tracks number of blocking groups
			lOuterindices = [] # Tracks atom indices part of blocking groups

			# Loop continues until an atom within the same grouping layer as atomSymbol is found to the left
			while LN not in self.__ATOMS or lBlockCounter > 0:
				if LN in self.__ATOMS:
					LNindex -= 1 # Decrement index because leftward expansion
				if LN == ')':
					lBlockCounter += 1
				if LN == '(' and lBlockCounter != 0: # LN can exist through parenthesis group (
					lBlockCounter -= 1
				if LN in self.__BONDS and lBlockCounter == 0: # Only add bond on same scope as atomHandler call
					LNbond = LN
				if LN in self.__NUMBERS and lBlockCounter == 0: # Only add number group on same scope as atomHandler call
					lOuterNumGroup = self.numbersHandler(LNPos)
				if LN in self.__ATOMS and lBlockCounter != 0:
					lOuterindices.append(LNindex)
				LNPos -= 1
				if LNPos < 0: # Breaks loop if scope is outside of SMILEScode range
					LNgroup = ""
					break
				LN = self.SMILES[LNPos]
			else:
				LNindex -= 1 # Decrement index by one to account for not evaluating final atom
				# Combine LN info and determine whichGroup
				if lOuterNumGroup: # If a number group existed, include in LNgroup
					LNindices = [LNindex, lOuterNumGroup[1], atomIndex]
					LNgroup = LN + '(' +  lOuterNumGroup[0] + ')' + LNbond + atomSymbol
				else: # If not, do not add any number group
					LNindices = [LNindex, atomIndex]
					LNgroup = LN + LNbond + atomSymbol
				# Determine the functional group information
				info = self.whichGroup(LNPos, atomPosition, LNgroup, LNindices)
				if info is not False:
					for group in info:
						FGinfo.append(group)

		# If there are entries in FGinfo, return it. Otherwise return False
		if FGinfo:
			return FGinfo
		else:
			return False # Returns False for no matches, otherwise the whichGroup info is returned based on atom expansion

	# Connects an explicit charge to the enclosed bracket group. Called when a + or - is encountered
	def chargeHandler(self,chargeSymbol, chargePostion, atomIndex):

		# Charges are simple within SMILEScode. Each individual charge is enclosed
		# in brackets with the charged atom, in the form [R(+,-)], like [C-] or [N+]. Charged hygroden
		# atoms are removed from the SMILEScode by the formatSMILEScode() function, since
		# no functioanl group needs to contain an expicitly charged hydrogen
		# SMILEScode cut from charge position to capture the form listed above,
		# and all functional groups containing that charge can be evaluted from it
		chargeGroup = self.SMILES[chargePostion-2:chargePostion+2]
		info = self.whichGroup(chargePostion-2, chargePostion+1, chargeGroup, [atomIndex])
		# If info is False, then False is returned. Likewise, groups will be returned if any are found
		return info # Returns False for no matches, otherwise the whichGroup info is returned based on a bracketed charge group

	# Called from within each handler to determine which functional groups the clutter of atoms could be apart of
	def whichGroup(self,startPosition, endPosition, group, atomindices):

		# Holds templates which contain group, sent to expandGroup to attempt a fullMatch
		portionMatches = []

		# Only one fullMatch group is possible
		fullMatch = False

		# Loop through group templates in FGlist
		for line in open('FGlist.txt', 'r'):

			# Variables
			lineInfo = re.compile(r'\S+').findall(line) # Grab line info seperated by whitespace
			FGtemplate = lineInfo[0].replace('[R]', 'R') # Replace with single R instead of bracketed R
			FGname = lineInfo[1] # Name comes after template on a single line
			difference = len(FGtemplate) - len(group) # Length difference to determine fullMatch or portion match

			match = self.checkGroup(group, FGtemplate, False) # Compare the group with the template, with partial matches allowed

			# If group is larger than the template, or no match was found, skip
			if difference < 0 or match is False:
				continue

			# Equivalent lengths with a True match means an identical match
			elif difference == 0 and match is True and fullMatch is False:
				fullMatch = True
				fullDict = self.createFGAtomDict(atomindices, FGtemplate)
				fullInfo = [FGname, FGtemplate, atomindices, fullDict]

			# Unequivalent lengths with a True match means a portion match
			elif difference > 0 and match is True:
				portionMatches.append([FGtemplate, FGname])

		# If a fullMatch was found without any portions, return group match information
		if fullMatch and not portionMatches:
			return [fullInfo]

		# If a fullMatch and some portions were found, determine if the portion groups can be completed with expandGroup.
		# Take expandGroup if it exists, otherwise take fullMatch
		elif fullMatch and portionMatches:
			expandedGroup = self.expandGroup(startPosition, endPosition, group, atomindices, portionMatches)
			if expandedGroup is not False:
				return expandedGroup
			else:
				# print("group was ", group, " starting at ", startPosition, " and ending at ", endPosition)
				# print("fullInfo ", fullInfo)
				return [fullInfo]

		# If there are only portions, determine if any of them can be completed with expand group
		elif fullMatch is False and portionMatches:
			expandedGroup = self.expandGroup(startPosition, endPosition, group, atomindices, portionMatches)
			if expandedGroup is not False:
				return expandedGroup
			else:
				return False

		# If there is no matches at all, return False
		elif fullMatch is False and not portionMatches:
			return False # Returns a list of all matched groups, i.e. [FGinfo1, FGinfo2... FGinfoN], where FGinfo itself is a list

	# Called from whichGroup to attempt to expand the clutter of atoms made from the hanlder into the functional groups determined by whichGroup 
	def expandGroup(self, startPosition, endPosition, group, atomindices, portionMatches, recursive=False):

		matches = [] # List which holds all matches, where each entry is an FGinfo list, i.e. [FGinfo1, FGinfo2... FGinfoN]
		portionCounter = -1 # Count the portion to obtain start position of group in template via GITposition
		for portionGroup in portionMatches:
			portionCounter += 1 # Increment counter

			# Initialize variables
			template = portionGroup[0] # FG template which is attempting to be expanded into
			FGname = portionGroup[1] # The name of the templated functioanl group
			numTemplateAtoms = len(self.__ATOMSREGEX.findall(template)) # Find number of atoms in template
			posInTemplate = self.determineShift(group, template) # Determine the start position of group in template
			requiredGroup = "" # String which holds the dynamically expanding group as the SMILEScode is evaluted with respect to the template and group variables
			leftRequired = template[0:posInTemplate] # String of characters requried to complete the fragmented group that exist to the left of group inside the template
			finalAtomindices = [] # List that holds the indicies of the atoms being expanded and appended to requiredGroup from inside the SMILEScode

			# Two expansion booleans that determine if the leftward/rightward expansion logics need to be executed
			# Unless group begins at 0 (leftExpand False) or ends at the length of template (rightExpand False), both expansion booleans will be True
			leftExpand = not self.checkGroup(group, template[0:len(group)], True)
			rightExpand = not self.checkGroup(group, template[posInTemplate:len(template)], True)

			# print("BEGIN")
			# print("\nIntial Condtions\n")
			# print("group = ", group)
			# print("template = ", template)
			# print("FGname = ", FGname)
			# print("atomindices = ", atomindices)
			# print("startPosition = ", startPosition)
			# print("endPosition = ", endPosition)
			# print("rightExpand = ", rightExpand)
			# print("leftExpand = ", leftExpand)

			# If leftExpand is not required, then group conditions already satisfied lefthand template. Set the group info to expandGroup variables
			if leftExpand is False:
				requiredGroup = group
				for index in atomindices:
					finalAtomindices.append(index)

			# Determine if this is a recursive call. If so then pop final index. See documentation for explanation
			numRequiredGroupAtoms = len(self.__ATOMSREGEX.findall(requiredGroup))
			if recursive is True and numRequiredGroupAtoms != len(atomindices):
				finalAtomindices.pop(-1)

			# LN expansion BEGIN

			# Varibles
			LNPos = startPosition
			LNIndex = atomindices[0]
			LNtempPos = posInTemplate

			# LN Loop

			# Loop continues until expansion is proven false or all lefthand required symbol from template are found in the SMILEScode
			# Required group is dynamically being created, meaning it is gaining more and more symbols as the loop continues
			# Since LNtempPos is a tracker of that dynamic position with respect to requiredGroup, an equality check between the two
			# to determine if the symbol was matched or not is possible. In other words, if the len(requiredGroup) == len(template[0:LNtempPos]),
			# Then the loop is sufficient. If however the lengths are unequal, then the position in the smiles code did not match the position
			# of the template. Therefore, the loop should cease. The loop continues until the len(leftRequired) is reached or until a match is found invalid
			while leftExpand is True:
				# Decrement counters for leftward expansion
				LNPos -= 1
				LNtempPos -= 1

				# Check if requirement was satisfied or SMILEScode is out of scope
				if len(requiredGroup) == len(leftRequired) or LNPos < 0:
					break

				# Set new SMILEScode and template symbol variables
				LN = self.SMILES[LNPos]
				LNtemp = template[LNtempPos]

				# six cases: LNtemp can be in ATOMS, R, PARENTHESIS, CHARGES, BONDS, BRACKETS
				# If LNtemp is in ATOMS, R, BONDS, BRACKET, and a PARENTHESIS or NUMBER is found, simply continue until the actual one is found
				if LNtemp == 'R':
					lBlockCounter = 0
					while LN not in self.__ATOMS or lBlockCounter > 0:
						if LN == ')':
							lBlockCounter += 1
						if LN == '(':
							lBlockCounter -= 1
						if lBlockCounter != 0 and LN in self.__ATOMS:
							LNIndex -= 1
						if (lBlockCounter == 0 and LN != "(" and LN not in self.__ATOMS and LN not in self.__NUMBERS) or LNPos < 0:
							break
						LNPos -= 1
						if LNPos < 0:
							break
						LN = self.SMILES[LNPos]
					else:
						LNIndex -= 1
						finalAtomindices.append(LNIndex)
						requiredGroup += LN

				elif LNtemp in self.__BRACKETS:
					# Many molecules do not have any charges, so check if there are any  in the smilescode before executing the loop
					# If there are none, the template will never match.
					numberOfBrackets = len(re.compile(r'\[').findall(self.SMILES))
					if numberOfBrackets > 0:
						lBlockCounter = 0
						while LN not in self.__BRACKETS or lBlockCounter > 0:
							if LN == ')':
								lBlockCounter += 1
							if LN == '(':
								lBlockCounter -= 1
							if lBlockCounter != 0 and LN in self.__ATOMS:
								LNIndex -= 1
							if (lBlockCounter == 0 and LN != '(' and LN not in self.__BRACKETS and LN not in self.__NUMBERS) or LNPos < 0:
								break
							LNPos -= 1
							if LNPos < 0:
								break
							LN = self.SMILES[LNPos]
						else:
							requiredGroup += LN

				elif LNtemp in self.__ATOMS:
					lBlockCounter = 0
					while LN in self.__NUMBERS or lBlockCounter > 0:
						if LN == ')':
							lBlockCounter += 1
						if LN == '(':
							lBlockCounter -= 1
						if lBlockCounter != 0 and LN in self.__ATOMS:
							LNIndex -= 1
						LNPos -= 1
						if LNPos < 0:
							break
						LN = self.SMILES[LNPos]
					else:
						if LN == LNtemp:
							LNIndex -= 1
							finalAtomindices.append(LNIndex)
							requiredGroup += LN

				elif LNtemp == ')':
					if LN in self.__NUMBERS:
						numGroupinfo = self.numbersHandler(LNPos)
						numGroup = '(' + numGroupinfo[0] + ')' # Leftward number is always outer group, never linear
						if numGroup == template[LNtempPos-len(numGroup)+1:LNtempPos+1]:
							requiredGroup += numGroup
					elif LN == ')':
						lBlockCounter = 1
						temp = LNtempPos
						tempSymbol = LNtemp
						requiredOuterGroup = ""
						while tempSymbol != '(':
							temp -= 1
							requiredOuterGroup += tempSymbol
							tempSymbol = template[temp]
						requiredOuterGroup = requiredOuterGroup[::-1]
						while lBlockCounter > 0:
							if LNPos in self.__ATOMS:
								LNindex -= 1
							if LN == ')':
								lBlockCounter += 1
							if LN == '(':
								lBlockCounter -= 1
							LNPos -= 1
							if LNPos < 0:
								break
							LN = self.SMILES[LNPos]
						else:
							outerGroup = self.SMILES[LNPos:LNPos+len(requiredGroup)-1] + ')'
							if outerGroup == requiredOuterGroup:
								requiredGroup += outerGroup

				elif LNtemp in self.__BONDS:
					if LN in self.__NUMBERS and (LNPos-1) >= 0:
						if self.SMILES[LNPos-1] == LNtemp:
							requiredGroup += self.SMILES[LNPos-1]
							LNPos -= 1
					if LN == LNtemp:
						requiredGroup += LN

				elif LNtemp in self.__CHARGES:
					if LN in self.__CHARGES:
						requiredGroup += LN

				elif LNtemp == '(':
					if LN == '(':
						requiredGroup += LN

				# If LNtemp was not satisfied by the SMILEScode, then the functional group does not exist in the SMILEScode
				# Check made at end of loop to ensure equality check was performed on the SMILEScode
				if len(requiredGroup) != len(template[LNtempPos:posInTemplate]):
					break

			# Continue to next possble functional group is loop failed
			if len(requiredGroup) != len(leftRequired) and leftExpand is True:
				# print("Left expansion failed for ", group, " into ", template)
				# print("END\n")
				continue
			# If there was a leftExpansion and the loop did not fail, then attatch the group and indices to finalAtomindices
			# Skipped if leftExpand was not required
			elif len(requiredGroup) == len(leftRequired) and leftExpand is True:
				finalAtomindices = finalAtomindices[::-1]
				requiredGroup = requiredGroup[::-1]
				requiredGroup += group
				for index in atomindices:
					finalAtomindices.append(index)

			# LN expansion END

			# RN expansion BEGIN

			# Variables
			RNPos = endPosition # RNPosition in SMILEScode
			RNIndex = atomindices[-1] # RNindex in SMILEScode
			RNtempPos = len(requiredGroup) - 1 # RNtemplateposition, intialized at final position of group within template
			outerExpandindices = [] # List of indices passed into recursive OuterExpansion call
			innerExpandindices = [] # List of indices passed into recursive InnerExpansion call
			InnerExpansion = OuterExpansion = 0 # Holds recursive information when called upon

			# If leftExpand was false, group length lends itself to RNtempPos
			if leftExpand is False:
				RNtempPos = len(group) - 1
			# print("\nIntial Right Expansion")
			# print("RNpos = ", RNPos)
			# print("RNtempPos = ", RNtempPos)
			# print("requiredGroup = ", requiredGroup)
			while rightExpand is True:
				# Increment counters for rightward expansion
				RNPos += 1
				RNtempPos += 1

				# If requiredGroup reaches the length of the template or the SMILEScode goes out of scope, break the loop
				if len(requiredGroup) == len(template) or RNPos > self.__SMILESlength -1 :
					break

				# Set new SMILEScode and template symbol variables
				RN = self.SMILES[RNPos]
				RNtemp = template[RNtempPos]
				# print("RNtemp is = ", RNtemp)
				# print("RN = ", RN)
				# print("requiredGroup = ", requiredGroup)
				# print("RNtempPos = ", RNtempPos)

				# If an outerGroup in the SMILEScode is encountered, collect all of its postional data to prepare for a recursive call
				if RN == '(':
					# Clear indices for a new outer group
					outerExpandindices.clear()
					innerExpandindices.clear()
					rBlockCounter = 1 # Intialized as 1 due to starting '('
					correlatingParenthPos = RNPos
					correlatingParenth = self.SMILES[RNPos]
					outerIndex = RNIndex
					innerInidicies = [outerIndex]
					while rBlockCounter > 0:
						correlatingParenthPos += 1
						correlatingParenth = self.SMILES[correlatingParenthPos]
						if correlatingParenth in self.__ATOMS:
							outerIndex += 1
							innerInidicies.append(outerIndex)
						if correlatingParenth == '(':
							rBlockCounter += 1
						if correlatingParenth == ')':
							rBlockCounter -= 1
					# Inner symbol is first symbol inside parenthesis, outerSymbol is first symbol after final parenthesis
					innerSymbol = self.SMILES[RNPos+1]
					innerIndex = RNIndex + 1
					outerSymbol = self.SMILES[correlatingParenthPos+1]
					outerIndex += 1
					for index in finalAtomindices:
						outerExpandindices.append(index)
						innerExpandindices.append(index)
					# Append the parenthesied group index i.e. C(R....)R, appending either the first R index in the parenthesis to expand outside, or vice versa
					if RNtemp == "(": # If RNtemp itself is an outer group, append the inner and outer indices as they will be required later on
						outerExpandindices.append(innerIndex)
						innerExpandindices.append(outerIndex)
					# Append the starting index at the end, see pop call
					outerExpandindices.append(outerIndex-1) # Ensure that OuterExpansion begins on the correct index, see pop() call at the top of expandGroup
					innerExpandindices.append(innerIndex-1) # Ensure that InnerExpansion begins on the correct index, see pop() call at the top of expandGroup

				if RNtemp == 'R':
					if RN == '(':
						# Only an R can exist after the parenthesis
						if innerSymbol in self.__ATOMS:
							InnerExpansion = self.expandGroup(LNPos, RNPos, requiredGroup, innerExpandindices, [portionGroup], True)
							if recursive is True:
								return InnerExpansion
							elif InnerExpansion is not False and recursive is False:
								matches.append(InnerExpansion)
						elif outerSymbol in self.__ATOMS:
							OuterExpansion = self.expandGroup(LNPos, correlatingParenthPos, requiredGroup, outerExpandindices, [portionGroup], True)
							if recursive is True:
								return OuterExpansion
							elif OuterExpansion is not False and recursive is False:
								matches.append(OuterExpansion)
						break
					if RN in self.__NUMBERS:
						numGroup = self.numbersHandler(RNPos)
						if numGroup[3] == "Linear":
							requiredGroup += numGroup[0]
							finalAtomindices.append(numGroup[1])
						elif (RNPos+1) < self.__SMILESlength:
							if self.SMILES[RNPos+1] in self.__ATOMS:
								# print(RNIndex)
								RNIndex += 1
								requiredGroup += self.SMILES[RNPos+1]
								RNPos += 1
								finalAtomindices.append(RNIndex)
					elif RN in self.__ATOMS:
						RNIndex += 1
						requiredGroup += RN
						finalAtomindices.append(RNIndex)


				elif RNtemp in self.__ATOMS:
					if RN == '(':
						# The specific atom can exist only after the parenthesis
						if innerSymbol == RNtemp:
							InnerExpansion = self.expandGroup(LNPos, RNPos, requiredGroup, innerExpandindices, [portionGroup], True)
							if recursive is True:
								return InnerExpansion
							elif InnerExpansion is not False and recursive is False:
								matches.append(InnerExpansion)
						elif outerSymbol == RNtemp:
							OuterExpansion = self.expandGroup(LNPos, correlatingParenthPos, requiredGroup, outerExpandindices, [portionGroup], True)
							if recursive is True:
								return OuterExpansion
							elif OuterExpansion is not False and recursive is False:
								matches.append(OuterExpansion)
						break
					if RN in self.__NUMBERS:
						numGroup = self.numbersHandler(RNPos)
						if numGroup[3] == "Linear" and numGroup[0] == RNtemp:
							requiredGroup += numGroup[0]
							finalAtomindices.append(numGroup[1])
						elif (RNPos+1) < self.__SMILESlength:
							if self.SMILES[RNPos+1] == RNtemp:
								RNIndex += 1
								requiredGroup += self.SMILES[RNPos+1]
								RNPos += 1
								finalAtomindices.append(RNIndex)
					elif RN == RNtemp:
						RNIndex += 1
						requiredGroup += RN
						finalAtomindices.append(RNIndex)

				elif RNtemp in self.__BRACKETS:
					if RN == '(':
						# Only an R can exist after the parenthesis
						if innerSymbol in self.__BRACKETS:
							InnerExpansion = self.expandGroup(LNPos, RNPos, requiredGroup, innerExpandindices, [portionGroup], True)
							if recursive is True:
								return InnerExpansion
							elif InnerExpansion is not False and recursive is False:
								matches.append(InnerExpansion)
						elif outerSymbol in self.__BRACKETS:
							OuterExpansion = self.expandGroup(LNPos, correlatingParenthPos, requiredGroup, outerExpandindices, [portionGroup], True)
							if recursive is True:
								return OuterExpansion
							elif OuterExpansion is not False and recursive is False:
								matches.append(OuterExpansion)
						break
					if RN in self.__NUMBERS:
						numGroup = self.numbersHandler(RNPos)
						if numGroup[3] == "Linear" and numGroup[0][0] == '[':
							loopCounter = 0
							for char in numGroup:
								loopCounter += 1
								if char != template[RNtempPos]:
									break
							if loopCounter == len(numGroup):
								requiredGroup += numGroup[0]
								finalAtomindices.append(numGroup[1])
								RNtempPos += (len(numGroup[0]) - 2)
						elif (RNPos+1) < self.__SMILESlength:
							if self.SMILES[RNPos+1] in self.__BRACKETS:
								requiredGroup += RN
								RNPos += 1
					elif RN in self.__BRACKETS:
						requiredGroup += RN

				elif RNtemp in self.__BONDS:
					if RN == '(':
						# Only an R can exist after the parenthesis
						if innerSymbol in self.__BONDS:
							InnerExpansion = self.expandGroup(LNPos, RNPos, requiredGroup, innerExpandindices, [portionGroup], True)
							if recursive is True:
								return InnerExpansion
							elif InnerExpansion is not False and recursive is False:
								matches.append(InnerExpansion)
						elif outerSymbol in self.__BONDS:
							OuterExpansion = self.expandGroup(LNPos, correlatingParenthPos, requiredGroup, outerExpandindices, [portionGroup], True)
							if recursive is True:
								return OuterExpansion
							elif OuterExpansion is not False and recursive is False:
								matches.append(OuterExpansion)
						break
					if RN in self.__NUMBERS and (RNPos + 1) < self.__SMILESlength:
						if self.SMILES[RNPos+1] == RNtemp:
							requiredGroup += self.SMILES[RNPos+1]
							RNPos += 1
					elif RN == RNtemp:
						requiredGroup += RN

				elif RNtemp in self.__CHARGES:
					if RN in self.__CHARGES:
						requiredGroup += RN

				elif RNtemp == '(':
					parenthesiedTemp = '('
					while RNtemp != ')':
						RNtempPos += 1
						RNtemp = template[RNtempPos]
						parenthesiedTemp += RNtemp
					numParenthesiedTempAtoms = len(self.__ATOMSREGEX.findall(parenthesiedTemp))
					# print("parenthesiedTemp = ", parenthesiedTemp)
					# print(numParenthesiedTempAtoms)
					if RN in self.__NUMBERS:
						numGroup = self.numbersHandler(RNPos)
						# If group is outer, add it. If it is linear, incompatable and group does not exist where number is represented at
						if numGroup[3] != "Outer":
							break
						if numParenthesiedTempAtoms == 1:
							outerParenthesied = '(' + numGroup[0] + ')'
							validNumGroup = self.checkGroup(outerParenthesied, parenthesiedTemp, True)
							# print("validNumGroup = ", validNumGroup)
							if validNumGroup is True:
								requiredGroup += outerParenthesied
								finalAtomindices.append(numGroup[1])
								# RNtempPos += (len(outerParenthesied) - 2s)
						elif numParenthesiedTempAtoms == 2 and numGroup[0] == parenthesiedTemp[1]: # If the atom in the number group is equal to the first atom inside the parenthesis, then continue
							validNumGroup = False
							# print("Hrere")
							LNrequired = parenthesiedTemp[2]
							lBlockCounter = 0
							LNIndex = numGroup[1]
							LNPos = numGroup[4] - 1
							LN = self.SMILES[LNPos]
							while LN != LNrequired or lBlockCounter > 0:
								if LN in self.__ATOMS and LNrequired == 'R' and lBlockCounter == 0:
									validNumGroup = True
									LNIndex -= 1
									break
								if LN == ')':
									lBlockCounter += 1
								if LN == '(':
									lBlockCounter -= 1
								if lBlockCounter != 0 and LN in self.__ATOMS:
									LNIndex -= 1
								if (lBlockCounter == 0 and LN not in self.__ATOMS) or LNPos < 0:
									validNumGroup = False
									break
								LNPos -= 1
								if LNPos < 0:
									validNumGroup = False
									break
								LN = self.SMILES[LNPos]
							else:
								validNumGroup = True
								LNIndex -= 1
							if validNumGroup is True:
								# print("Found valid numGroup")
								outerParenthesied = '(' + numGroup[0] + LN + ')'
								requiredGroup += outerParenthesied
								finalAtomindices.append(numGroup[1])
								finalAtomindices.append(LNIndex)
						elif numParenthesiedTempAtoms == 0:
							print("Parenthesied temp has 0 atoms, fatal error.")
							print(self.SMILES)
							sys.exit()
					elif RN == '(':
						innerParenthesied = self.SMILES[RNPos:RNPos+len(parenthesiedTemp)-1] + ')'
						outerParenthesied = '(' + self.SMILES[correlatingParenthPos+1:correlatingParenthPos+len(parenthesiedTemp)-1] + ')'
						innerExpand = self.checkGroup(innerParenthesied, parenthesiedTemp, True)
						outerExpand = self.checkGroup(outerParenthesied, parenthesiedTemp, True)
						# print("innerParenthesied = ", innerParenthesied)
						# print("outerParenthesied = ", outerParenthesied)
						if innerExpand is True:
							numAtoms = len(self.__ATOMSREGEX.findall(innerParenthesied))
							if numAtoms == 2: # If there was 2 atoms, append another index right before the last
								outerExpandindices.insert(len(outerExpandindices)-1,innerIndex+1)
								# print("Outer expanding with ", requiredGroup+innerParenthesied, " and the index's", indices)
							OuterExpansion = self.expandGroup(startPosition, correlatingParenthPos, requiredGroup+innerParenthesied, outerExpandindices, [portionGroup], True)
							if recursive is True:
								return OuterExpansion
							elif OuterExpansion is not False and recursive is False:
								matches.append(OuterExpansion)
						elif outerExpand is True:
							numAtoms = len(self.__ATOMSREGEX.findall(outerParenthesied))
							if numAtoms == 2: # If there was 2 atoms, append another index right before the last
								innerExpandindices.insert(len(innerExpandindices)-1,outerIndex+1)
							# print("Inner expanding with ", requiredGroup+ outerParenthesied, " and the indices ", indices)
							InnerExpansion = self.expandGroup(startPosition, RNPos, requiredGroup+outerParenthesied, innerExpandindices, [portionGroup], True)
							if recursive is True:
								return InnerExpansion
							elif InnerExpansion is not False and recursive is False:
								matches.append(InnerExpansion)
						break

				elif RNtemp == ')':
					if RN in self.__NUMBERS and (RNPos+1) < self.__SMILESlength:
						if self.SMILES[RNPos+1] == ')':
							requiredGroup += self.SMILES[RNPos+1]
							RNPos+=1
					elif RN == ')':
						requiredGroup += RN

				# If the symbol(s) were not added, then the functional group does not match in the SMILEScode. Break
				if len(requiredGroup) != len(template[0:RNtempPos+1]):
					# print(requiredGroup, " != ", template[0:RNtempPos+1])
					# print(RNtemp, " did not match ", RN)
					break
				# else:
				# 	print("\nSymbols matched ")
				# 	print(requiredGroup, " = ", template[0:RNtempPos+1])
				# 	print("\n")

			# If the expansion failed, or if a recursive call was successfull, continue to the next group
			if OuterExpansion is 0 and InnerExpansion is 0: # If rightExpand stopped on a symbol by symbol analysis at the highest stack call
				# print("Determining the outcome of the expansion logic")
				finalCheck = self.checkGroup(requiredGroup, template, True)
				if finalCheck is True:
					# print("finalCheck = ", True)
					# print(requiredGroup, " = ", template[0:RNtempPos+1])
					finalDict = self.createFGAtomDict(finalAtomindices, template)
					# print("finalCheck = ", finalCheck)
					# print("recursive = ", recursive)
					# print("len(portionMatches) = ", len(portionMatches))
					FGInfo = [FGname, template, finalAtomindices, finalDict]
					if recursive is True and len(portionMatches) == 1:
						return FGInfo
					else:
						matches.append(FGInfo)
				if finalCheck is False:
					# print("finalCheck = ", False)
					# print(requiredGroup, " != ", template[0:RNtempPos+1])
					if recursive is True:
						return False
					else:
						continue
			elif recursive is False: # If a decision was made by expansion logics, no need to check anything. Simply move to the next grou
				# print("A decision was made on Expansion logics")
				# print("InnerExpansion = ", InnerExpansion)
				# print("OuterExpansion = ", OuterExpansion)
				continue
			else:
				print("Fatal error ")
				print("finalCheck = ", finalCheck)
				print("finalAtomindices = ", finalAtomindices)
				print("numTemplateAtoms = ", numTemplateAtoms)
				print("group = ", group)
				print("requiredGroup = ", requiredGroup)
				print("template = ", template)
				self.printLists()

		# print("\nEXPAND GROUP RESULTS FOR THE GROUP ", group)
		if matches:
			# for group in matches:
				# print(group)
			return matches
		else:
			# print("FOUND NOTHING!")
			return False

	# Connects the two atoms related by number, or ring. Determines the atom and bond associated with that number
	def numbersHandler(self, numPosition):

		for RINGGROUP in self.RINGPOSITIONS:

			# If the number OPENS the group, return the CLOSING info
			if RINGGROUP[0][0] == numPosition:
				return(RINGGROUP[1][1], RINGGROUP[1][2], RINGGROUP[1][3], RINGGROUP[0][4], RINGGROUP[1][0])

			# If the number CLOSES the group, return the OPENING info
			elif RINGGROUP[1][0] == numPosition:
				return(RINGGROUP[0][1], RINGGROUP[0][2], RINGGROUP[0][3], RINGGROUP[1][4], RINGGROUP[0][0])

	# Creates a dictionary of index to symbol from the template and indicies of a functional group
	def createFGAtomDict(self, FGindices, FGtemplate):
		FGdict = []
		indexCounter = -1
		#print(FGindices)
		for symbol in FGtemplate:
			if symbol in self.__ATOMS or symbol == 'R':
				indexCounter += 1
				#print(indexCounter)
				FGdict.append([symbol, FGindices[indexCounter]])
		#print(FGdict)
		return FGdict

	# Checks if a group is equal to template, or if a group is a portion of a template, Boolean
	def checkGroup(self, group, template, full=True):

		group = group.upper()
		difference = len(template) - len(group)
		# ##print("Checking ", group, " aginast ", template)
		# ##print("There is a difference of ", difference, " while full = " , full)
		# Group larger than template, return False
		if difference < 0:
			return False

		# If only fullMatch is possible, but difference exists, return False
		if difference > 0 and full is True:
			return False

		# Group same length as template
		if difference == 0:
			groupCounter = -1
			for symbol in template:
				groupCounter += 1
				groupSymbol = group[groupCounter]
				if groupSymbol not in self.__ATOMS and symbol == 'R': # If symbol is an atom where it should be
					return False
				elif groupSymbol != symbol and symbol != 'R': # If symbol is not an R-group, then it must be an identical component in template
					return False
			return True

		# Group smaller than template
		if difference > 0 and full is False:
			for shift in range(0, difference+1):
				templateCounter = - 1 + shift
				symbolCounter = - 1
				for symbol in group:
					symbolCounter += 1
					templateCounter += 1
					templateSymbol = template[templateCounter]
					if templateSymbol != symbol and templateSymbol != 'R':
						break
					if symbol not in self.__ATOMS and templateSymbol == 'R':
						break
					if symbolCounter == len(group) - 1:
						return True
			return False

	# Determines where in the template a given group appears. Called within expandGroup only, gauranteed to return a value
	def determineShift(self,group, template):
		group = group.upper()
		difference = len(template) - len(group)
		for shift in range(0, difference+1):
			templateCounter = - 1 + shift
			symbolCounter = - 1
			for symbol in group:
				symbolCounter += 1
				templateCounter += 1
				templateSymbol = template[templateCounter]
				if templateSymbol != symbol and templateSymbol != 'R':
					break
				if symbol not in self.__ATOMS and templateSymbol == 'R':
					break
				if symbolCounter == len(group) - 1:
					return shift
