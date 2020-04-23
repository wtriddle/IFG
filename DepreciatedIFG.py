# Written by William Riddle from 9/23/19 - 4/1/20

import re
import sys
from dataManipulators import dataManipulators

class ifg(dataManipulators):
	
	def __init__(self, inputSMILES):

		self.__ALCOHOLINFO = []
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
		self.CHARGEDMOLECULE = True if len(re.compile(r'\+\-').findall(inputSMILES)) != 0 else False

		self.SMILES = self.formatSMILEScode(inputSMILES)
		self.__SMILESlength = len(self.SMILES)
		self.__AROMATICCOUNT = 0
		self.__NONAROMATICCOUNT = 0
		self.__RINGCOUNT = 0
		self.AROMATICINDICES = self.initializeAROMATICINDICES()
		self.CYCLICINDICES = self.initializeCYCLICINDICES()
		self.RINGPOSITIONS = self.initializeRINGPOSITIONS()
		self.determineRingCounts()
		self.SMILES = self.SMILES.upper()
		self.FGDATAFINAL = self.identifyFunctionalGroups()

	def identifyFunctionalGroups(self):

		# Collection Phase of script is run first to determine all possible functional groups, placed into __FGdata
		self.collectFunctionalGroups()

		print(self.__FGdata)

		# Evaluation Phase to evalute collected groups, finalize them with accuracy and precision, and gather other various data parameters
		data = self.evaluateFunctionalGroups()

		return data
	
	def collectFunctionalGroups(self):
			
		atomIndex = -1

		for SMILEScodePos, symbol in enumerate(self.SMILES):
			
			if symbol in self.__ATOMS:
				atomIndex += 1

				directGroup = self.whichGroup(symbol, SMILEScodePos, atomIndex)
				if directGroup:
					for group in directGroup:
						self.__FGdata.append(group)
				del(directGroup) 

				if symbol == 'O':
					# Determine if an alochol is present at a given oxygen
					self.determineAlcoholGroup(SMILEScodePos, atomIndex)
				elif symbol == 'N':
					# Determine if a primaryAmine is present at a given nitrogen
					self.determinePrimaryAmine(SMILEScodePos, atomIndex)

	def evaluateFunctionalGroups(self):
			
		# Scrub the FGdata list for removal of incorrect groups or repeated groups
		temp1 = self.repetitionScrub()
		temp2 = self.alcoholScrub(temp1)
		temp3 = self.heirarchyScrub(temp2)
		temp4 = self.amideScrub(temp3)

		for index,group in enumerate(self.__ALCOHOLINFO):
			temp4.append(['Alcohol', 'ROH', [group[0],self.__ALCOHOLICINDICES[index]]])

		allFGstemp = []
		for innerGroup in temp4:
			allFGstemp.append(innerGroup[:])

		acetals = self.determineAcetalGroups(allFGstemp)
		for acetal in acetals:
			allFGstemp.append(acetal[:])
			temp4.append(acetal[:])

		# print(temp4)

		FGdataFinaltemp = self.containmentScrub(temp4)

		alcoholCount = len(self.__ALCOHOLINFO)

		FGdataFinal = self.determineCyclicGroups(FGdataFinaltemp)
		allFGs = self.determineCyclicGroups(allFGstemp)

		print("FGdataFinal = ", FGdataFinal)
		print("allFGs = ", allFGs)

		FGdataFinalDict = self.createFGDataDict(FGdataFinal)
		allFGsDict = self.createFGDataDict(allFGs)
		FGdataFinalDict.update({"totalAlcohols" : alcoholCount})
		allFGsDict.update({"totalAlcohols" : alcoholCount})

		# Occurance of [NH+],[NH2+], or [NH3+] gaurnetees amino acid strucutre
		if self.AMINOACID is True:
			FGdataFinalDict.update({'AminoAcid' : 1})
		else:
			FGdataFinalDict.update({'AminoAcid' : 0})

		if self.__NONAROMATICCOUNT + self.__AROMATICCOUNT == self.__RINGCOUNT:
			FGdataFinalDict.update({'RingCount' : self.__RINGCOUNT})
			FGdataFinalDict.update({'aromaticRingCount' : self.__AROMATICCOUNT})
			FGdataFinalDict.update({'nonAromaticRingCount' : self.__NONAROMATICCOUNT})
		else:
			print("Fatal ring count error. Exiting..")
			self.printLists()
			print(self.__NONAROMATICCOUNT, " + ", self.__AROMATICCOUNT, " != ", self.__RINGCOUNT )
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
					self.__ALCOHOLICINDICES.append(atomIndex)
					self.__ALCOHOLINFO.append([LNIndex-1])
					self.__ALCOHOLINFO[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
					self.__ALCOHOLINFO[-1].append("FinalO")

		# First atom alcohol case
		# Check for oxygen at first position of SMILEScode, and a single carbon bound to it
		elif atomIndex == 0 and self.SMILES[SMILEScodePos+1].upper() == 'C':
			self.__ALCOHOLICINDICES.append(atomIndex)
			self.__ALCOHOLINFO.append([atomIndex+1]) # Automatically append next index if satisfied
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
			self.__ALCOHOLINFO[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
			self.__ALCOHOLINFO[-1].append("First O")
		# (O) single oxygen outer groups are automatically alcoholic oxygens since they are isolated
		elif self.SMILES[SMILEScodePos-1:SMILEScodePos+2] == '(O)':
			self.__ALCOHOLICINDICES.append(atomIndex)
			self.__ALCOHOLINFO.append([atomIndex-1]) # Atom index pointing to R group bound to (O), i.e. C(O), pointing to the C
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
			self.__ALCOHOLINFO[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
			self.__ALCOHOLINFO[-1].append("Beginning (")
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
				self.__ALCOHOLICINDICES.append(atomIndex)
				self.__ALCOHOLINFO.append([LNIndex-1])
				self.__ALCOHOLINFO[-1].append(numRGroups) # Append the number of Rgroups to newly appended alcohol group
				self.__ALCOHOLINFO[-1].append("Closing )")
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
			for alocholIndex, acetalCarbon in enumerate(self.__ALCOHOLINFO):
				if acetalCarbon[0] in groupIndices: # Check if the acetal carbon is part of the ether. If so, it is part of an acetal group
					numRgroups = acetalCarbon[1] # Grab the numRGroups from associated list
					if numRgroups == 3: # 3 Rgropus means Hemiacetal group
						# Create Hemiacetal group from known information
						HemiacetalIndices = []
						for index in groupIndices: # All ether indices are part of Hemiacetal indices
							HemiacetalIndices.append(index)
						HemiacetalIndices.append(self.__ALCOHOLICINDICES[alocholIndex]) # As well as alcoholic oxygen
						distinctAcetal = not all(index in acetalGroupindices for index in HemiacetalIndices) # Check if the Hemiacetal group has already been found
						if distinctAcetal is True: # If not already found, then it is a new Hemiacetal group. Add its contents to acetalGroups
							if HemiacetalIndices not in acetalGroupindices:
								acetalGroups.append(['Hemiacetal', '', HemiacetalIndices])
								for index in HemiacetalIndices:
									acetalGroupindices.append(index)
								acetalGroupindices.append(HemiacetalIndices)
						break
					elif numRgroups == 4: # 4 Rgroups means Hemiketal group
						# Create Hemiketal group from known information
						HemiketalIndices = []
						for index in groupIndices: # All ether indices are part of Hemiketal indices
							HemiketalIndices.append(index)
						HemiketalIndices.append(self.__ALCOHOLICINDICES[alocholIndex]) # As well as alcoholic oxygen
						distinctAcetal = not all(index in acetalGroupindices for index in HemiketalIndices) # Check if the Hemiketal group has already been found
						if distinctAcetal is True: # If not already found, then it is a new Hemiketal group. Add its contents to acetalGroups
							if HemiketalIndices not in acetalGroupindices:
								acetalGroups.append(['Hemiketal', '', HemiketalIndices])
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
					# Add a distinctly found acetal group to indices
					distinctAcetal = not all(index in acetalGroupindices for index in acetalIndices) # Check if the acetal group has already been found
					if distinctAcetal is True:
						acetal = ["Acetal", "", acetalIndices]
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

						# Append the opening and closing (or correlating) information into a single list.
						# Each list within RINGPOSITIONS indicates one set of ring information
						openingInfo = [SMILEScodePos, openingAtom, atomIndex, openingAtomPos]
						closingInfo = [correlatingPos, correlatingAtom, correlatingIndex, correlatingAtomPos]
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
		reFormatted = ""
		for SMILEScodePos, symbol in enumerate(SMILEScode):
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

		atomIndex = -1
		for SMILEScodePos, symbol in enumerate(self.SMILES):

			if symbol in self.__ATOMS:
				atomIndex += 1

			# If a new ring has been found in the SMILES code via a number 
			if symbol in self.__NUMBERS and SMILEScodePos not in evaluatedNumbers:
				# print("Found ", symbol, " at position, ", SMILEScodePos, " with an index of ", atomIndex )
				# print("\n")

				scopeIndices = [[]] # Index by the scope of paretnehsis groups. I.e. 0th scope is scope where number was found, 1st is an inner parethensis, etc.
				rBlockCounter = 0
				evaluatedNumbers.append(SMILEScodePos)
				scopeIndices[0].append(atomIndex)
				RNPos = SMILEScodePos + 1
				RNindex = atomIndex
				RN = self.SMILES[RNPos]

				# Symbol is equal to number which opened the ring. Loop untl the same number is encountered by RN
				while RN != symbol:

					if RN in self.__ATOMS:
						RNindex+=1
						if len(scopeIndices) == rBlockCounter: # If a new nested parenthesis or deeper parenthesis is found, i.e. in between a parenthesis group (...()..,)
							scopeIndices.append([RNindex]) # Create a new list tracking the inner indices of that nested parenthesis group
						else: # If the parenthesis scope is not deeper, then access the current scope and append the atom index to it
							scopeIndices[rBlockCounter].append(RNindex)

					if RN == '(':
						rBlockCounter += 1
					if RN == ')':
						del(scopeIndices[rBlockCounter])
						rBlockCounter -= 1
					RNPos+=1
					RN = self.SMILES[RNPos]

				else:
					evaluatedNumbers.append(RNPos) 
					for scope in scopeIndices:
						for index in scope: # Add all indices within that scope to CYCLICindices
							if index not in CYCLICINDICES:
								CYCLICINDICES.append(index)

		return CYCLICINDICES
		
	# Create the AROMATICINDICES list
	def initializeAROMATICINDICES(self):
		AROMATICINDICES = []
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

		formattedData = dataList[:]
		groupCounter = - 1

		while groupCounter <= len(formattedData):
			groupCounter += 1 

			if groupCounter == len(formattedData):
				break

			group = formattedData[groupCounter]
			groupIndices = group[2]

			# If the group contians an alcohol
			if len(re.compile(r'acid|oxide|oxime').findall(group[0])) != 0:
				# And it does not contain an alcoholic oxygen index
				if not any(index in self.__ALCOHOLICINDICES for index in groupIndices):
					# Remove its contents
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


	def amideScrub(self, dataList):
		
		formattedData = dataList[:]

		groupCounter = -1
		while groupCounter <= len(formattedData):
			
			groupCounter+=1
			if groupCounter == len(formattedData):
				break

			group = formattedData[groupCounter]
			groupName = group[0]
			groupIndices = group[2]

			for compareGroupCounter, compareGroup in enumerate(formattedData):
				
				if compareGroupCounter == groupCounter:
					continue
					
				compareName = compareGroup[0]	
				compareIndices = compareGroup[2]

				if compareName == groupName == "Amide":
					
					if all(index in groupIndices for index in compareIndices):
						print("Popping ", compareGroup, " from ", group)
						formattedData.pop(compareGroupCounter)
						groupCounter = -1
						break

					elif all(index in compareIndices for index in groupIndices):
						print("Popping ", group, " from ", compareGroup)
						formattedData.pop(groupCounter)
						groupCounter = -1
						break
		return formattedData

	# Scrubs FG dataList for full containment of groups. i.e. Ester takes precedence over ketones and ethers
	def containmentScrub(self, dataList):

		# Initializations
		formattedData = dataList[:] # Copy input dataList
		print("FormattedData = ", formattedData)
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
			groupName = group[0]
			groupTemplate = group[1]
			groupIndices = group[2]
			compareGroupCounter = -1 # Set compareGroupCounter for each new group


			# print("Determing Containment for ", group)
			# Compare this group with all other groups found in SMILES for containment
			for compareGroup in formattedData:

				compareGroupCounter += 1 # Track index
				compareGroupName = compareGroup[0]
				compareTemplate = compareGroup[1]
				compareIndices = compareGroup[2]

				# Determine if all of the group atoms are contained within the compareGroup
				# i.e. check for full group containment, like a Ketone contained within an Ester
				# All of its indices must be within a larger groups indices
				fullGroupContainment = all(index in compareIndices for index in groupIndices)

				# If group is fully contained within compareGroup, that is every atom is in comparegroup, remove the group and restart loop without that group
				if fullGroupContainment and len(groupIndices) < len(compareIndices):
					# print("Popping G", group, " from ", compareGroup)
					formattedData.pop(groupCounter)
					groupCounter = -1
					break

				# Similalry, check if comapreGroup satisfies this condition
				fullCompareGroupContainment = all(index in groupIndices for index in compareIndices)

				# If compareGroup is fully contained within a group, that is every atom is in group, remove the group and restart loop without that group
				if fullCompareGroupContainment and len(compareIndices) < len(groupIndices):
					# print("Popping CG", compareGroup, " from ", group)
					formattedData.pop(compareGroupCounter)
					groupCounter = -1
					break

				# Same atoms, but one more specific than other based on R groups
				numMainGroupAtoms = len(self.__ATOMSREGEX.findall(groupTemplate)) - len(self.__RGROUPREGEX.findall(groupTemplate))
				numMainCompareGroupAtoms = len(self.__ATOMSREGEX.findall(compareTemplate)) - len(self.__RGROUPREGEX.findall(compareTemplate))

				if fullGroupContainment and fullCompareGroupContainment:
					
					if numMainCompareGroupAtoms > numMainGroupAtoms:
						# print("Popping G", group, " from ", compareGroup)
						formattedData.pop(groupCounter)
						groupCounter = -1
						break

					elif numMainCompareGroupAtoms < numMainGroupAtoms:
						# print("Popping CG", compareGroup, " from ", group)
						formattedData.pop(compareGroupCounter)
						groupCounter = -1
						break

		return formattedData 

	# Determines if two functional groups are exactly identical
	def isGroupInList(self, dataList, group):

		if dataList:
			for compareGroup in dataList:

				compareTemplate = compareGroup[1]
				compareIndices = compareGroup[2]

				groupTemplate = group[1]
				groupIndices = group[2]

				# Indices must be equal
				if all(index in groupIndices for index in compareIndices):
				
					groupDict = self.createFGAtomDict(groupIndices,groupTemplate)
					compareDict = self.createFGAtomDict(compareIndices,compareTemplate)

					index = 0
					for atom in groupDict:
						if atom not in compareDict:
							break
						index+=1
					
					if index == len(groupDict) == len(compareDict):
						return True
					else: 
						continue
		return False 

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
		print("ALCOHOLICINDICES = ", self.__ALCOHOLINFO)
		print("CYCLICINDICES = ", self.CYCLICINDICES)
		print("AROMATICINDICES = ", self.AROMATICINDICES)
		print("RINGPOSITIONS = ", self.RINGPOSITIONS)
	
	# Called from within each handler to determine which functional groups the clutter of atoms could be apart of
	def whichGroup(self, symbol, symbolPosition, atomIndex):

		matches = []
		# Loop through group templates in FGlist
		for line in open('FGlist.txt', 'r'):

			# Variables
			lineInfo = re.compile(r'\S+').findall(line) # Information seperated by a white space
			lineInfo[0] = lineInfo[0].replace('[R]','R')
			FGtemplate = lineInfo[0]
			# Must be a chrage in SMILES for charged FG to match
			if len(re.compile(r'\+\-').findall(FGtemplate)) != 0 and self.CHARGEDMOLECULE is False:
				continue
			
			if symbol in FGtemplate:
				
				expansionPoints = [] # String positions from which the template can be expanded into the SMILES code

				for templatePosition, templateSymbol in enumerate(FGtemplate):
					if templateSymbol == symbol:
						expansionPoints.append(templatePosition)

				for expansionPoint in expansionPoints:
					expandedGroup = self.expandGroup(symbolPosition, symbolPosition, symbol, [atomIndex], FGtemplate, expansionPoint, False)
					if expandedGroup != False:
						expandedGroup.insert(0, lineInfo[1]) # Insert FG name 
						matches.append(expandedGroup)

				del(expandedGroup)

		if matches:
			return matches
		else:
			return False

	# Called from whichGroup to attempt to expand the clutter of atoms made from the hanlder into the functional groups determined by whichGroup 
	def expandGroup(self, startPosition, endPosition, group, atomindices, template, expansionPoint, recursive=False):

		# Initialize variables
		if len(group) == 1:
			posInTemplate = expansionPoint
		else:
			posInTemplate = self.determineShift(group, template)
		requiredGroup = "" # Tracking the SMILEScode template equivalent.
		leftRequiredTemplate = template[0:posInTemplate] 
		finalAtomindices = [] 

		leftExpand = not recursive
		# if template == "RN(=O)=O":
		# 	print("BEGIN")
		# 	print("\nIntial Condtions\n")
		# 	print("group = ", group)
		# 	print("template = ", template)
		# 	print("FGname = ", FGname)
		# 	print("atomindices = ", atomindices)
		# 	print("startPosition = ", startPosition)
		# 	print("endPosition = ", endPosition)
		# 	print("leftExpand = ", leftExpand)
		# 	print("posInTemplate = ", posInTemplate)

		if leftExpand is False:
			requiredGroup = group
			for index in atomindices:
				finalAtomindices.append(index)

		# Left Expansion
		LNPos = startPosition
		LNIndex = atomindices[0]
		LNtempPos = posInTemplate
		# if template == "RN(=O)=O":
		# 	print("\nLeft Expand Inital")
		# 	print("LNPos = ", LNPos)
		# 	print("LNIndex = ", LNIndex)
		# 	print("LNtempPos = ", LNtempPos)
		


		while leftExpand is True:
			
			LNPos -= 1
			LNtempPos -= 1

			# Check if requirement was satisfied or SMILEScode is out of scope
			if self.checkGroup(requiredGroup, leftRequiredTemplate) or LNPos < 0:
				break

			LN = self.SMILES[LNPos]
			LNtemp = template[LNtempPos]

			# Maintin LN expansion on the same scope from where it was called from
			# Sets position of LN exactly one position to the left of '(' when the next iteration of the loop is incremented on
			if LN == ')':
				
				lBlockCounter = 1

				while lBlockCounter > 0 or LN == ')':
					LNPos -= 1
					LN = self.SMILES[LNPos]
					if LN in self.__ATOMS:
						LNIndex -= 1
					elif LN == '(':
						lBlockCounter -= 1
					elif LN == ')':
						lBlockCounter += 1
				else:
					# Do not skip template symbol. Add one to keep constant when next iteration is increment on
					LNtempPos += 1
					continue

			# Pass invidiual opening parenthesis if LNtemp is not inside parenthesis group
			elif LN == '(' and LNtemp not in ['(',')']:
				LNtempPos+=1 
				continue

			# Template switch case analysis on SMILEScode
			if LNtemp == 'R':
				
				# Find a non numerical symbol to the left
				while LN in self.__NUMBERS:
					LNPos -= 1
					if LNPos < 0:
						break
					LN = self.SMILES[LNPos]

				# Check if that symbol is an atom
				else: 
					if LN in self.__ATOMS:
						LNIndex -= 1
						finalAtomindices.insert(0,LNIndex)
						requiredGroup = LN + requiredGroup

			elif LNtemp in self.__BRACKETS:
				
				# Find a non numerical symbol to the left
				while LN in self.__NUMBERS:
					LNPos -= 1
					if LNPos < 0:
						break
					LN = self.SMILES[LNPos]

				# Check if that symbol is a bracket
				else:
					if LN in self.__BRACKETS:
						requiredGroup = LN + requiredGroup

			elif LNtemp in self.__ATOMS:
				
				# Find non numerical symbol to the left
				while LN in self.__NUMBERS:
					LNPos -= 1
					if LNPos < 0:
						break
					LN = self.SMILES[LNPos]
				
				# Check if that symbol is identical to the template symbol 
				else:
					if LN == LNtemp:
						LNIndex -= 1
						finalAtomindices.insert(0,LNIndex)
						requiredGroup = LN + requiredGroup
			
			elif LNtemp in self.__BONDS:
							
				# Find a non numerical symbol to the left
				while LN in self.__NUMBERS:
					LNPos -= 1
					if LNPos < 0:
						break
					LN = self.SMILES[LNPos]

				# Check if that symbol is a bond
				else:
					if LN in self.__BONDS:
						requiredGroup = LN + requiredGroup


			elif LNtemp == ')':
				
				# Retrieve parenthesied group within template
				holder = LNtempPos
				while template[LNtempPos] != '(':
					LNtempPos -= 1
				parenthesiedTemp = template[LNtempPos:holder+1]

				# Numerical symbol in SMILES can satisfy parenthesied groups in template
				if LN in self.__NUMBERS:
					numGroupinfo = self.numbersHandler(LNPos)
					numGroup = '(' + numGroupinfo[0] + ')' # Leftward number is always outer group, never linear
					if self.checkGroup(numGroup, parenthesiedTemp):
						requiredGroup = numGroup + requiredGroup
						finalAtomindices.insert(0,numGroupinfo[1])
						continue

				# SMILES parenthesis groups cut into by the appropriate length can satisfy parenthesied groups in template
				if self.SMILES[LNPos+1] == '(':
					parenthesiedSMILES = self.SMILES[LNPos+1:LNPos+len(parenthesiedTemp)] + ')'
					if self.checkGroup(parenthesiedSMILES, parenthesiedTemp):
						requiredGroup = parenthesiedSMILES + requiredGroup
						if len(self.__ATOMSREGEX.findall(parenthesiedSMILES)) == 2:
							finalAtomindices.insert(0,LNIndex)
							finalAtomindices.insert(1,LNIndex+1)
						else:
							finalAtomindices.insert(0,LNIndex)
						LNPos+=1

			elif LNtemp in self.__CHARGES and LN in self.__CHARGES:
					requiredGroup = LN + requiredGroup

			elif LNtemp == '(' and LN == '(':
					requiredGroup = LN + requiredGroup

			# If the symbol required in the template did not match the symbol in the SMILEScode
			if self.checkGroup(requiredGroup, template[LNtempPos:posInTemplate]) is False:
				return False
		
		if leftExpand:
			requiredGroup += group
			finalAtomindices.append(atomindices[0])


		# Right Expansion
		RNPos = endPosition 
		RNIndex = atomindices[-1] 
		if recursive:
			RNIndex = expansionPoint
		RNtempPos = len(requiredGroup) - 1 # Intialize at final position of group within template

		# if template == "RN(=O)=O":
		# 	print("\nRight Expand Initial:")
		# 	print("RNPos = ", RNPos)
		# 	print("RNIndex = ", RNIndex)
		# 	print("RNtempPos = ", RNtempPos)

		while self.checkGroup(requiredGroup,template[0:RNtempPos+1]):
			
			RNPos += 1
			RNtempPos += 1

			# If requiredGroup reaches the length of the template or the SMILEScode goes out of scope, break the loop
			if len(requiredGroup) == len(template) or RNPos > self.__SMILESlength -1 :
				break

			RN = self.SMILES[RNPos]
			RNtemp = template[RNtempPos]

			# if template == "RN(=O)=O":
			# 	print("\nRight Expand Conditons:")
			# 	print("RNPos = ", RNPos)
			# 	print("RNIndex = ", RNIndex)
			# 	print("RNtempPos = ", RNtempPos)
			# 	print("RN = ", RN)
			# 	print("RNtemp = ", RNtemp)
			# 	print("requiredGroup = ", requiredGroup)

			if RN == '(':

				outerExpandindices = finalAtomindices[:]
				innerExpandindices = finalAtomindices[:] 
				rBlockCounter = 1 
				openParPos = closeParPos = RNPos # openParenthesisPosition, closeParenthesisPosition
				innerIndexInit = RNIndex # Points one index to left of first index inside of parentehsis group
				closePar = self.SMILES[RNPos]

				# Loop ceases on final SMILES parenthesis
				while rBlockCounter > 0:
					closeParPos += 1
					closePar = self.SMILES[closeParPos]
					if closePar in self.__ATOMS:
						RNIndex += 1
					if closePar == '(':
						rBlockCounter += 1
					if closePar == ')':
						rBlockCounter -= 1
					
				outerIndexInit = RNIndex # Points to final index in parenthesis group

				# Inner symbol is first symbol inside SMILES parenthesis, outerSymbol is first symbol after final SMILES parenthesis
				innerSymbol = self.SMILES[openParPos+1]
				outerSymbol = self.SMILES[closeParPos+1]

				if RNtemp in self.__ATOMS or RNtemp in self.__BONDS or RNtemp in self.__BRACKETS:

					if innerSymbol == RNtemp:
						InnerExpansion = self.expandGroup(LNPos, openParPos, requiredGroup, innerExpandindices, template, expansionPoint=innerIndexInit, recursive=True)
						if InnerExpansion:
							return InnerExpansion

					if outerSymbol == RNtemp:
						OuterExpansion = self.expandGroup(LNPos, closeParPos, requiredGroup, outerExpandindices, template, expansionPoint=outerIndexInit, recursive=True)
						if OuterExpansion:
							return OuterExpansion

				elif RNtemp == 'R':

					if innerSymbol in self.__ATOMS:
						InnerExpansion = self.expandGroup(LNPos, openParPos, requiredGroup, innerExpandindices, template, expansionPoint=innerIndexInit, recursive=True)
						if InnerExpansion:
							return InnerExpansion

					if outerSymbol in self.__ATOMS:
						OuterExpansion = self.expandGroup(LNPos, closeParPos, requiredGroup, outerExpandindices, template, expansionPoint=outerIndexInit, recursive=True)
						if OuterExpansion:
							return OuterExpansion

				elif RNtemp == '(':
					
					# Retrieve parenthesied group within template
					holder = RNtempPos
					while template[RNtempPos] != ')':
						RNtempPos += 1
					parenthesiedTemp = template[holder:RNtempPos+1]

					# Determine expansion paths based on SMILES and templated parenthesis groups
					innerParenthesied = self.SMILES[openParPos:openParPos+len(parenthesiedTemp)-1] + ')'
					outerParenthesied = '(' + self.SMILES[closeParPos+1:closeParPos+len(parenthesiedTemp)-1] + ')'
					innerExpand = self.checkGroup(outerParenthesied, parenthesiedTemp)
					outerExpand = self.checkGroup(innerParenthesied, parenthesiedTemp)

					if outerExpand is True:
						numAtoms = len(self.__ATOMSREGEX.findall(innerParenthesied))
						for innerIndex in range(1,numAtoms+1):
							outerExpandindices.append(innerIndexInit+innerIndex)
						OuterExpansion = self.expandGroup(startPosition, closeParPos, requiredGroup+innerParenthesied, outerExpandindices, template, expansionPoint=outerIndexInit, recursive=True)
						if OuterExpansion:
							return OuterExpansion

					if innerExpand is True:
						numAtoms = len(self.__ATOMSREGEX.findall(outerParenthesied))
						for outerIndex in range(1,numAtoms+1):
							innerExpandindices.append(outerIndexInit+outerIndex)
						InnerExpansion = self.expandGroup(startPosition, openParPos, requiredGroup+outerParenthesied, innerExpandindices, template, expansionPoint=innerIndexInit, recursive=True)
						if InnerExpansion:
							return InnerExpansion
				
				return False

			if RNtemp == 'R':
				if RN in self.__NUMBERS:
					numGroup = self.numbersHandler(RNPos)
					if numGroup[0] in self.__ATOMS and RNtempPos == len(template) - 1:
						requiredGroup += numGroup[0]
						finalAtomindices.append(numGroup[1])
					elif (RNPos+1) < self.__SMILESlength:
						if self.SMILES[RNPos+1] in self.__ATOMS:
							RNIndex += 1
							requiredGroup += self.SMILES[RNPos+1]
							RNPos += 1
							finalAtomindices.append(RNIndex)
				elif RN in self.__ATOMS:
					RNIndex += 1
					requiredGroup += RN
					finalAtomindices.append(RNIndex)


			elif RNtemp in self.__ATOMS:
				if RN in self.__NUMBERS:
					numGroup = self.numbersHandler(RNPos)
					if RNtempPos == len(template) - 1 and numGroup[0] == RNtemp:
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
				if RN in self.__NUMBERS:
					numGroup = self.numbersHandler(RNPos)
					if numGroup[0][0] == '[':
						if self.SMILES[RNPos:RNPos+len(numGroup)] == numGroup[0]:
							requiredGroup += numGroup[0]
							finalAtomindices.append(numGroup[1])
							RNtempPos += (len(numGroup[0]) - 1)
					elif (RNPos+1) < self.__SMILESlength:
						if self.SMILES[RNPos+1] in self.__BRACKETS:
							requiredGroup += RN
							RNPos += 1
				elif RN in self.__BRACKETS:
					requiredGroup += RN

			elif RNtemp in self.__BONDS:
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
				
				holder = RNtempPos				
				while template[RNtempPos] != ')':
					RNtempPos += 1
				parenthesiedTemp = template[holder:RNtempPos+1]
				numParenthesiedTempAtoms = len(self.__ATOMSREGEX.findall(parenthesiedTemp))

				if RN in self.__NUMBERS:
					
					numGroup = self.numbersHandler(RNPos)

					if numParenthesiedTempAtoms == 1:
						numParenthesied = '(' + numGroup[0] + ')'
						if self.checkGroup(numParenthesied, parenthesiedTemp):
							requiredGroup += numParenthesied
							finalAtomindices.append(numGroup[1])

					elif numParenthesiedTempAtoms == 2 and numGroup[0] == parenthesiedTemp[1]:
						
						LNrequired = parenthesiedTemp[2]
						lBlockCounter = 0
						LNIndex = numGroup[1]
						LNPos = numGroup[3] - 1
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
			elif RNtemp == ')':
				if RN in self.__NUMBERS and (RNPos+1) < self.__SMILESlength:
					if self.SMILES[RNPos+1] == ')':
						requiredGroup += self.SMILES[RNPos+1]
						RNPos+=1
				elif RN == ')':
					requiredGroup += RN


		finalCheck = self.checkGroup(requiredGroup, template)
		if finalCheck:
			FGInfo = [template, finalAtomindices]
			# if template == "RN(=O)=O":
			# 	print(FGInfo)
			return FGInfo
		else:
			return False

	# Connects the two atoms related by number, or ring. Determines the atom and bond associated with that number
	def numbersHandler(self, numPosition):

		# List in form [atomSymbol, atomIndex, atomPosition, correlatingNumberPosition]
		for ring in self.RINGPOSITIONS:

			# If the number OPENS the group, return the CLOSING info
			if ring[0][0] == numPosition:
				return(ring[1][1], ring[1][2], ring[1][3], ring[1][0])

			# If the number CLOSES the group, return the OPENING info
			elif ring[1][0] == numPosition:
				return(ring[0][1], ring[0][2], ring[0][3], ring[0][0])

	# Checks if a group is equal to template
	def checkGroup(self, group, template):

		group = group.upper()
		difference = len(template) - len(group)

		if difference != 0:
			return False

		elif difference == 0:
			groupCounter = -1
			for symbol in template:
				groupCounter += 1
				groupSymbol = group[groupCounter]
				if groupSymbol not in self.__ATOMS and symbol == 'R': 
					return False
				elif groupSymbol != symbol and symbol != 'R': 
					return False
			return True

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

	# Finds the number of aromatic, nonaromatic, and total ring counts
	def determineRingCounts(self):

		self.__RINGCOUNT = len(self.RINGPOSITIONS)

		# Simplification of aromatic/nonaromatic count
		allAtoms = ''.join(self.__ATOMSREGEX.findall(self.SMILES))

		if allAtoms.islower():
			self.__AROMATICCOUNT = len(self.RINGPOSITIONS)
			return 0

		elif allAtoms.isupper():
			self.__NONAROMATICCOUNT = len(self.RINGPOSITIONS)
			return 0
		
		# Upper is nonaromatic, lower is aromatic
		for ring in self.RINGPOSITIONS:
			print("Evaluating ring ", ring)
			
			# Both are aromatic
			if ring[0][1].islower() and ring[1][1].islower():

				openingNumberPos = ring[0][0]

				if self.SMILES[openingNumberPos+1] in self.__ATOMS:
					
					if self.SMILES[openingNumberPos+1].islower():
						self.__AROMATICCOUNT+=1

					elif self.SMILES[openingNumberPos+1].isupper():
						self.__NONAROMATICCOUNT+=1

				else:
					self.__AROMATICCOUNT+=1

			# Both are non aromatic
			elif ring[0][1].isupper() and ring[1][1].isupper():
				
				openingNumberPos = ring[0][0]

				if self.SMILES[openingNumberPos+1] in self.__ATOMS:
					
					if self.SMILES[openingNumberPos+1].islower():
						self.__AROMATICCOUNT+=1

					elif self.SMILES[openingNumberPos+1].isupper():
						self.__NONAROMATICCOUNT+=1

				else:
					self.__NONAROMATICCOUNT+=1
					
			# They are differnt, must be nonaromatic
			else:
				print("Found non aromatic ring ")
				self.__NONAROMATICCOUNT+=1
		
	def createFGAtomDict(self, FGindices, FGtemplate):
		FGdict = []
		indexCounter = -1
		for symbol in FGtemplate:
			if symbol in self.__ATOMS or symbol == 'R':
				indexCounter += 1
				FGdict.append([symbol, FGindices[indexCounter]])
		return FGdict
