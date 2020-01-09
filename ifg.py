# Written by William Riddle from 9/23/19 -

import re
import sys

# Notes about variables and SMILEScode
# ATOMS and RGROUP are identical, but they are different for clarity of program when determining groups.
# The RGROUPREGEX is a regex just for the letter R, as it appeares in the FGlist templates.

# Globals
RINGPOSITIONS = []
FGdata = []
GITposition = []
CYCLICINDICES = [] # Tracks all cyclic indicies
AROMATICINDICIES = [] # Tracks only aromatic cyclic indicies

# Constants
ATOMS = ['C','O','N','c','n','o']
RGROUP = ['C','O','N','c','n','o']
CHARGES = ['+','-']
PARENTHESIS = ['(', ')']
BRACKETS = ['[', ']']
BONDS = ['=', '#']
NUMBERS = ['1','2','3','4','5','6','7','8','9']

# Reuglar Expressions
ATOMSREGEX = re.compile(r'[a-zA-Z]')
NUMBERSREGEX = re.compile(r'\d')
RGROUPREGEX = re.compile(r'R')

def ifg(SMILES, SMILESWH):

	# Initializations
	global SMILEScode
	SMILEScode = formatSMILEScode(SMILES) # Definition of SMILEScode
	global SMILEScodelength
	SMILEScodelength = len(SMILEScode) # Definition of SMILEScodelength
	SMILEScodePos = atomIndex = -1 # Counters
	initializeRINGPOSITIONS() # Create dynamic model of SMILEScode
	initializeCYCLICINDICES() # Find the cyclic indices of the SMILEScode

	# Main Loop
	for symbol in SMILEScode:

		SMILEScodePos += 1

		# GITposition is symbol specific, reset its contents each new symbol
		GITposition.clear()

		if symbol in ATOMS:
			atomIndex += 1 # Track the index

		if symbol in PARENTHESIS:
			# Do nothing, skip to next symbol
			continue

		elif symbol in CHARGES:
			# #print("SMILEScode[", SMILEScodePos, "] = ", symbol, " and belongs to the CHARGE group")
			# Call chargeHandler on the charge and pass in positional parameters
			chargeGroup = chargeHandler(symbol, SMILEScodePos, atomIndex)
			# print("chargeHandler found ", chargeGroup, "\n")
			if chargeGroup is not False:
				for group in chargeGroup:
					FGdata.append(group)
			del(chargeGroup)

		# Any atom except for the first and last atom may be sent to atomHandler
		elif symbol in ATOMS and SMILEScodePos != 0 and SMILEScodePos != SMILEScodelength - 1:
			#print("SMILEScode[", SMILEScodePos, "] = ", symbol, " and belongs to the ATOMS group, and is atom number ", atomIndex)
			# Call atomHandler on the atom, pass in positional parameters
			directGroup = atomHandler(symbol, SMILEScodePos, atomIndex)
			# print("atomHandler found ", directGroup, "\n")
			if directGroup is not False:
				# Direct group may have multiple FG's found, independently add them all
				for group in directGroup:
					FGdata.append(group)
			del(directGroup)

		elif symbol in BONDS:
			##print("SMILEScode[", SMILEScodePos, "] = ", symbol, " and belongs to the BONDS group")
			# Call bondHandler on the bond, pass in positional parameters
			bondGroup = bondHandler(symbol, SMILEScodePos, atomIndex)
			# print("bondHandler found ", bondGroup, "\n")
			if bondGroup is not False:
				for group in bondGroup:
					FGdata.append(group)
			del(bondGroup)

	# for group in FGdata:
	# 	print(group, " \n")
	FGdataFinal = evaluateFGdata() # Validate and format FGdata to capture all present groups
	determineCyclicGroups(FGdataFinal) # Find the cyclic/aromatic functionalGroups if any exist
	FGdataDict = {}

	# Loop through FGdataFinal and increment the determined functioanl groups to find their number of occurances
	for group in FGdataFinal:
		keys = []
		for key in FGdataDict.keys():
			keys.append(key)
		if group[0] not in keys:
			FGdataDict.update({group[0] : 1})
		else:
			FGdataDict[group[0]] += 1
		del(keys) # Delete the data after usage to prevent potential data leaks

	# Use hydrogen-explict SMILEScode to determine alcohol count
	alcoholCount = len(re.compile(r'OH').findall(SMILESWH))
	if alcoholCount != 0:
		FGdataDict.update({'Alcohol' : alcoholCount})

	# Occurance of [NH+],[NH2+], or [NH3+] gaurnetees amino acid strucutre
	# print(len(re.compile(r'\[[nN]H[23]?\+\]').findall(SMILESWH)))
	AminoAcid = True if len(re.compile(r'\[[nN]H[23]?\+\]').findall(SMILESWH)) != 0 else False
	if AminoAcid is True:
		FGdataDict.update({'AminoAcid' : 1})
	else:
		FGdataDict.update({'AminoAcid' : 0})

	# Determine the number of aromatic and nonAromatic rings from RINGPOSITIONS global
	ringCount = len(RINGPOSITIONS)
	aromaticRingCount = nonAromaticRingCount = 0
	# print("RINGPOSITIONS = ", RINGPOSITIONS)
	for ring in RINGPOSITIONS:
		if ring[0][1].islower(): # Lowercase means aromatic ring in SMILEScode
			aromaticRingCount+=1
		if ring[0][1].isupper(): # Uppercase means non-aromatic ring in SMILEScode
			nonAromaticRingCount+=1

	# Make sure the total ring count is valid, then add the data to FGdataDict
	if nonAromaticRingCount + aromaticRingCount == ringCount:
		FGdataDict.update({'RingCount' : ringCount})
		FGdataDict.update({'aromaticRingCount' : aromaticRingCount})
		FGdataDict.update({'nonAromaticRingCount' : nonAromaticRingCount})
	else:
		print("Fatal ring count error. Exiting..")
		print(nonAromaticRingCount, " + ", aromaticRingCount, " != ", ringCount )
		print(RINGPOSITIONS)
		print(SMILEScode)

	# Dump the data to prevent any data leaks or potential errors
	# print(RINGPOSITIONS)
	GITposition.clear()
	FGdata.clear()
	RINGPOSITIONS.clear()
	CYCLICINDICES.clear()
	AROMATICINDICIES.clear()

	# Return the functional group data for the given SMILEScode
	return (FGdataFinal, FGdataDict)


def initializeRINGPOSITIONS():

	# Initializations
	SMILEScodePos = atomIndex = -1 # Counters
	evaluatedNumbers = [] # List that tracks the string positions of the numbers whose ring junction information has already been determined

	# Main Loop
	for symbol in SMILEScode:

		# Two properoties which distinguish if a number is an outer group, i.e. ([R]), or a linear group, i.e. [R].
		# See documentation for explanation of this SMILEScode analysis
		openingProperty = ""
		closingProperty = ""
		##print(RINGPOSITIONS)
		SMILEScodePos += 1

		# Keep track of ATOM positional data for when a number is found
		if symbol in ATOMS:
			atomIndex += 1
			# Atom found before number is "opening" atom of the ring.
			openingAtom = symbol
			openingAtomPos = SMILEScodePos

		# When an unevalutaed number, or essentially a new ring, is encountered, evalute it
		if symbol in NUMBERS and SMILEScodePos not in evaluatedNumbers:
			chargeGroup = "" # Holds charged atoms if a chargeGroup opens or closes a ring. Reset to "" for subsequent charge groups
			##print("Found ", symbol, " at position ", SMILEScodePos)
			##print(evaluatedNumbers)

			# If the symbol before the number is not an atom, then it must be a bracket surrounding a charged atom.
			# Otherwise, it is always an atom and this special charge case will be skipped

			# The openingAtomPos variable points to the charged atom in the following case
			# The following if statement below is a special case where charged atom is opening the ring.
			# A closing bracket must precede the position of the number.
			if SMILEScode[SMILEScodePos-1] == ']' and SMILEScodePos >= openingAtomPos + 2:
				specialCounter = SMILEScodePos - 1 # Initialized at position of left bracket
				specialSymbol = SMILEScode[specialCounter]
				while specialSymbol != '[':
					chargeGroup += specialSymbol
					specialCounter -= 1
					# If the counter goes out of the scope of the SMILEScode, exit the script
					if specialCounter == -1:
						print("Fatal error in the creation of RINGPOSITIONS opening charge group on the smilescode ", SMILEScode)
						sys.exit()
					specialSymbol = SMILEScode[specialCounter]
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
			for subSymbol in SMILEScode[SMILEScodePos+1:SMILEScodelength]:
				correlatingPos += 1

				# Keep track of ATOM positional data for when the correlating atom is found
				if subSymbol in ATOMS:
					correlatingAtom = subSymbol
					correlatingIndex += 1
					correlatingAtomPos = correlatingPos

				# If the same number which began the sub loop is found, then the ring closure has been found
				# The correlatingAtom positional data and opening atom data will be paired together in RINGPOSITIONS
				if subSymbol == symbol:
					##print("Correlating number ", subSymbol, " was found at position ", correlatingPos)

					# If the ring opens with a charge, it must close with a charge. Re-assign positional variables
					# This special case is evaluated in the same manner as in the opening case
					if chargeGroup != "":
						correlatingChargeGroup = ""
						specialCounter = correlatingPos -1
						specialSymbol = SMILEScode[specialCounter]
						while specialSymbol != '[':
							correlatingChargeGroup += specialSymbol
							specialCounter -= 1
							# If the counter goes out of the scope of the SMILEScode, exit the script
							if specialCounter == -1:
								print("Fatal error in the creation of RINGPOSITIONS correlating charge group on the smilescode ", SMILEScode)
								sys.exit()
							specialSymbol = SMILEScode[specialCounter]
						correlatingChargeGroup = correlatingChargeGroup[::-1] # Flip because of leftward expansion

						# Charge group must be a length of at most 4. Anything larger suggests an error, so take the other atom first
						if len(correlatingChargeGroup) <= 4:
							# Reset the variables in this special case
							correlatingAtom = correlatingChargeGroup
							correlatingAtomPos = specialCounter

					# Determine linear or outer distinctions for each atom

					# Opening property distinction
					RNPos = openingAtomPos
					RN = SMILEScode[RNPos]
					while openingProperty == "":
						RNPos+=1
						if RNPos >= SMILEScodelength - 1:
							openingProperty = "Linear"
							break
						RN = SMILEScode[RNPos]
						if RN in ATOMS or RN in BONDS:
							openingProperty = "Outer"
							break
						if RN == ')':
							openingProperty = "Linear"
							break

					# Correlating or closing property distinction
					RNPos = correlatingPos
					RN = SMILEScode[RNPos]
					while closingProperty == "":
						RNPos+=1
						if RNPos >= (SMILEScodelength - 1):
							closingProperty = "Linear"
							break
						RN = SMILEScode[RNPos]
						if RN in ATOMS or RN in BONDS:
							closingProperty = "Outer"
							break
						if RN == ')':
							closingProperty = "Linear"
							break

					RINGPOSITIONS.append([
						[
						SMILEScodePos,
						openingAtom,
						atomIndex,
						openingAtomPos,
						openingProperty],
						[correlatingPos,
						correlatingAtom,
						correlatingIndex,
						correlatingAtomPos,
						closingProperty]
						])

					# Add numbers to evaluted list
					evaluatedNumbers.append(SMILEScodePos)
					evaluatedNumbers.append(correlatingPos)
					break

				# Error caught if a correlatingPos runs out of the scope of the SMILEScode
				if correlatingPos == SMILEScodelength - 1:
					print("There has been an error in the creation of the RINGPOSITIONS global \n")
					print("The SMILEScode code ", SMILEScode, " had an error")
					sys.exit()
	return 0


def evaluateFGdata():

	# Initializations
	formattedFGdata = [] # Captures the correct groups from FGdata
	groupCounter = -1 # Counter
	# print("FG data is ", FGdata)

	# Intial loop removes reprititons of certain FG's in FGdata
	# The result is placed into formattedFGdata
	for group in FGdata:
		# print(group)
		inFormattedData = isGroupInList(formattedFGdata, group)
		if inFormattedData is True:
			continue
		else:
			formattedFGdata.append(group)

	# Continue until formattedFGdata is looped through without any revisions
	# Removal of groups causes a reset of the loop.
	while groupCounter <= len(formattedFGdata):

		# Initializations
		groupCounter += 1 # Counter
		# If the loop limit has been reached, then exit the loop
		if groupCounter == len(formattedFGdata):
			break
		# Grab the group from groupCounter and set list contents to explicit varibales
		group = formattedFGdata[groupCounter]
		groupTemplate = group[1]
		groupindices = group[2]
		groupDict = group[3]
		compareGroupCounter = -1 # Set compareGroupCounter for each new group

		# Determine number of Maingroup atoms in group
		numMainGroupAtoms = len(groupindices) - len(RGROUPREGEX.findall(groupTemplate))

		# compareGroup Loop
		for compareGroup in formattedFGdata:
			compareGroupCounter += 1 # Counter

			# Grab the compareGroup from compareGroupCounter and set list contents to explicit varibales
			compareTemplate = compareGroup[1]
			compareindices = compareGroup[2]
			compareDict = compareGroup[3]

			# Determine number of Maingroup atoms in compareGroup
			numMainCompareAtoms = len(compareindices) - len(RGROUPREGEX.findall(compareTemplate))

			# Determine if all of the group atoms are contained within the compareGroup
			# i.e. check preemptively for group containment, for example a Ketone contained within an Ester
			fullGroupContainment = all(index in compareindices for index in groupindices)

			# If group is fully contained within a compareGroup, remove the group and restart loop without that group
			if fullGroupContainment is True and len(groupindices) < len(compareindices):
				formattedFGdata.pop(groupCounter)
				groupCounter = -1
				break

			# Similalry, check if comapreGroup satisfies this condition
			fullCompareGroupContainment = all(index in groupindices for index in compareindices)

			# If group is fully contained within a compareGroup, remove the group and restart loop without that group
			if fullCompareGroupContainment is True and len(compareindices) < len(groupindices):
				formattedFGdata.pop(compareGroupCounter)
				groupCounter = -1
				break

			# Determine mainGroupCrossover if there is no full containments
			mainGroupCrossover = 0
			for atom in groupDict:
			 	for compareAtom in compareDict:
					# If the atoms are symbol equivalent and both are main group atoms
			 		if atom[0] == compareAtom[0] and atom[0] != 'R' and compareAtom[0] != 'R':
						# If the atoms are the same index within the SMILEScode
			 			if int(atom[1]) == int(compareAtom[1]):
							# Then atoms must be of the exact same maingroup in the structue
			 				mainGroupCrossover += 1
			# If there is no relation, then group is independent of compareGroup, even if there is Rgroup crossover
			if mainGroupCrossover == 0:
				continue

			# If there is full mainGroupCrossover between the two groups, i.e. they have identical mainGroupAtoms
			if mainGroupCrossover == numMainCompareAtoms == numMainGroupAtoms:

				# If compareindices has more R groups, take it over group
				if len(groupindices) < len(compareindices):
					formattedFGdata.pop(groupCounter)
					groupCounter = -1
					break

				# If groupindices has more R groups, take it over compareGroup
				if len(groupindices) > len(compareindices):
					formattedFGdata.pop(compareGroupCounter)
					groupCounter = -1
					break

			# If all main group atoms in compareGroup are contained within group and both groups are the same length, but group has more main group atoms, then remove compareGroup
			if mainGroupCrossover == numMainCompareAtoms and numMainGroupAtoms > numMainCompareAtoms:
				formattedFGdata.pop(compareGroupCounter)
				groupCounter = -1
				break

			# If all main group atoms in group are contained within compareGroup and both groups are the same length, but compareGroup has more main group atoms, then remove group
			if mainGroupCrossover == numMainGroupAtoms and numMainCompareAtoms > numMainGroupAtoms:
				formattedFGdata.pop(groupCounter)
				groupCounter = -1
				break

	# Returned as FGdataFinal in IFG head function
	return formattedFGdata


def isGroupInList(formattedFGdata, group):

	# If there is an atom in formattedFGdata
	if formattedFGdata:
		# ##print("formattedFGdata is ", formattedFGdata)

		for compareGroup in formattedFGdata:
			# print("Compare group is ", compareGroup)
			# print("group is ", group)
			# Comparative Variables
			compareGroupTemplate = compareGroup[1]
			compareGroupindices = compareGroup[2]
			compareGroupDict = compareGroup[3]
			numMainCompareAtoms = len(compareGroupindices) - len(RGROUPREGEX.findall(compareGroupTemplate))

			# Group Variables
			groupTemplate = group[1]
			# print(SMILEScode)
			groupindices = group[2]
			groupDict = group[3]
			numMainGroupAtoms = len(groupindices) - len(RGROUPREGEX.findall(groupTemplate))

			# Check for equivalent number of atoms
			if len(compareGroupindices) == len(groupindices):
				mainGroupCrossover = 0
				# print(groupDict, groupTemplate)
				# print(compareGroupDict, compareGroupTemplate)
				for atom in groupDict:
					for compareAtom in compareGroupDict:

						# Check if atoms are maingroup equivalent, and if they are in the same positional index.
						# The atoms in both functionalGroup scopes are "identical" if so, so add to mainGroupCrossover
						if atom[0] == compareAtom[0] and atom[0] != 'R':
							if int(atom[1]) == int(compareAtom[1]):
								mainGroupCrossover += 1

				# If there is full maingroups crossover and the templates are the same length
				if mainGroupCrossover == numMainCompareAtoms == numMainGroupAtoms and len(compareGroupTemplate) == len(groupTemplate):
						# Then the exact group is already contained within formattedFGdata
						return True
			else: # if not equivalent lengths, continue
				continue
	return False


def bondHandler(bondSymbol, bondPostion, atomIndex):

	# Initialize positional variables to create group
	LNPos = RNPos = bondPostion # Initialized at the poistion of the bond
	LNindex = RNindex = atomIndex # Initialized at the position of the previous atom to the left of where the bond appears
	LN = RN = bondSymbol # Initialized as bondSymbol
	#print("LN = ", LN)


	# Special case which occurs commonly as R(=R)
	if SMILEScode[LNPos-1] == '(' and SMILEScode[LNPos-2] in ATOMS and SMILEScode[RNPos+1] in ATOMS and SMILEScode[RNPos+2] == ')':
		atomindices = [atomIndex, atomIndex+1]
		bondGroup = SMILEScode[LNPos-2] + "(=" + SMILEScode[RNPos+1] + ')'
		# ##print("Bond group from special case is ", bondGroup)
		FGinfo = whichGroup(LNPos-2, RNPos+2, bondGroup, atomindices)
		if FGinfo is not False:
			return FGinfo
		else:
			return False


	# LN loop
	lBlockCounter = 0
	lOuterParenthGroup = ""
	lOuterNumGroup = ""
	lOuterindices = []
	LNindex+=1 # Increment index to allow loop to account for offset of intialized atomIndex, which is already the lefthand index.
	# Therefore, must increment so that when an atom is encounter in the loop is actually points to the correct atom

	# Find the first non blocking group atom connected to the bond
	# If LN is an atom, loop not run
	while LN not in ATOMS and lBlockCounter > 0:
		LNPos -= 1
		LN = SMILEScode[LNPos]
		if LN in CHARGES and lBlockCounter == 0:
			return False
		if LN == ')':
			lBlockCounter += 1
		if LN == '(' and lBlockCounter != 0:
			lBlockCounter -= 1
		if LN in NUMBERS and lBlockCounter == 0:
			lOuterNumGroup = numbersHandler(LNPos)
		if lBlockCounter != 0:
			lOuterParenthGroup += LN
		if LN in ATOMS:
			LNindex -= 1
		if lBlockCounter != 0 and LN in ATOMS:
			lOuterindices.append(LNindex)

	# RN loop. If RN is already an atom, loop not run
	# The atom next to the bond symbol is always the atom involved in the bond, no blocking groups occur
	while RN not in ATOMS:
		RNPos += 1
		RN = SMILEScode[RNPos]
		if RN in BRACKETS:
			return False
		if RN in ATOMS:
			RNindex += 1

	# lOuterParenthGroup evaluation to see if the group should be included in bondGroup
	lOuterParenthGroup = lOuterParenthGroup[::-1] # Flip because appended from the left, so group appears backwards
	atomCount = len(ATOMSREGEX.findall(lOuterParenthGroup))
	numCount = len(NUMBERSREGEX.findall(lOuterParenthGroup))
	if atomCount > 1 or numCount != 0:
		lOuterParenthGroup = ""

	# Bond Group creation according to the specific case
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
	#print("Looking for group = ", bondGroup)
	FGinfo = whichGroup(LNPos, RNPos, bondGroup, atomindices)
	return FGinfo


def atomHandler(atomSymbol, atomPosition, atomIndex):
	# ##print("atomHandler is creating initial functional group shells...\n")
	# Initialize Variables
	LNPos = atomPosition - 1
	RNPos = atomPosition + 1
	RN = SMILEScode[RNPos]
	LN = SMILEScode[LNPos]
	# If atom handler called from within a charged and bracketed symbol, return False
	if RN in CHARGES:
		return False
	RNindex = LNindex = atomIndex
	RNbond = LNbond = ""
	FGinfo = []

# !!!!!!!!!!!!!RN LOOP BEGIN!!!!!!!!!!!!!!!!!

	# Methodology for loop in atomHandler differers from bondHandler and expandGroup.
	# Symbols and positions must be initalized next to the atom of interest because
	# the loop cannot begin on the atomSymbol itself, otherwise the loops would
	# never collect the necessary group information. Therefore, the RNPos variable
	# is incremented at the bottom of the loop instead of at the top and the information
	# is analyzed in the next run of the loop as opposed to the immediate loop run
	# ##print("-----RN LOOP-----")
	# ##print("INTIAL CONDITIONS:\n")
	# ##print("rBlockCounter = ", rBlockCounter)
	# ##print("rInnerParenthGroup = ", rInnerParenthGroup, "\n")
	# ##print("SMILEScode[", RNPos, "] = RN = ", RN)
	# ##print("RNindex = ", atomIndex)
	rBlockCounter = 0
	rOuterGroups = 0
	rInnerParenthGroup = "("
	rInnerParenthindices = []
	rInnerParenthGroupPositions = []
	numGroup = ""

	# Loop continues until an atom within the same grouping layer as atomSymbol is found to the right
	loopCounter = 0
	rOuterBoolean = True
	while rBlockCounter > 0 or RN not in ATOMS:
		##print(rOuterGroups)
		if RNPos == SMILEScodelength - 1:
			RN = ""
			##print("NO Right group exists...")
			break
		loopCounter += 1
		if RN in ATOMS:
			RNindex += 1
		# Add RN numbers handler
		# ##print(loopCounter)
		# IF RN in the same grouping layer as atomSymbol

		# If scope of RN loop is within a parenthesis, capture the information inside of it with respect to atomSymbol
		if rBlockCounter != 0:
			rInnerParenthGroup += RN
			rInnerParenthGroupPositions.append(RNPos)
			if rOuterBoolean is True:
				rOuterGroups += 1
				rOuterBoolean = False
			if RN in ATOMS:
				rInnerParenthindices.append(RNindex)
			if RN in NUMBERS:
				numGroup = numbersHandler(RNPos)
				rInnerParenthindices.append(numGroup[1])
				if numGroup[3] == "Linear":
					rInnerParenthGroup += numGroup[0]
				elif numGroup[3] == "Outer":
					outerGroup = '(' + numGroup[0] + ')'
					rInnerParenthGroup += outerGroup
		if rOuterBoolean is False and rBlockCounter == 0:
			##print("Reset rOuterBoolean to ", rOuterBoolean)
			rOuterBoolean = True
			##print("Found inner group and added one to ", rOuterGroups, " and set rOuterBoolean to ", rOuterBoolean)

		# ##print("RN LOOP INFO:\n")
		# ##print("rBlockCounter = ", rBlockCounter)
		# ##print("RNpos = ", RNPos)
		# ##print("SMILEScode[", RNPos, "] = RN = ", RN)
		# ##print("rInnerParenthGroup = ", rInnerParenthGroup)
		# ##print("RNindex = ", RNindex)
		# ##print("numGroup = ", numGroup)
		if RN in BONDS and rBlockCounter == 0:
			RNbond = RN
		if RN == '(':
			rBlockCounter += 1
		if RN == ')':
			rBlockCounter -= 1

		# If RN on the edge of a parenthesis, no RN exists, cease the loop
		if rBlockCounter < 0:
			RN = ""
			break

		RNPos += 1
		RN = SMILEScode[RNPos]
		# ##print("Made it to the end of this loop!")
		# ##print("RN = ", RN)
		# ##print("InnerLoop = ", innerLoop)
	else:
		RNindex += 1
		##print("CEASED ON ", loopCounter)
		##print("\nRN LOOP RESULT:\n")
		##print("rBlockCounter = ", rBlockCounter)
		##print("RNpos = ", RNPos)
		##print("SMILEScode[", RNPos, "] = RN = ", RN)
		##print("rInnerParenthGroup = ", rInnerParenthGroup)
		##print("RNindex = ", RNindex)
		##print("numGroup = ", numGroup)
		##print("Not terminated by a break!")

# !!!!!!!!!!!!! RN LOOP END !!!!!!!!!!!!!!!!!!!

	rInnerAtoms = ATOMSREGEX.findall(rInnerParenthGroup)

	# ##print("Matches = ", matches)
	atomCount = len(rInnerAtoms)

	##print("atomCount == ", atomCount)
	##print(rOuterGroups)

	# If there are more than 2 atoms inside the inner RN group and there is only one parenthesis group
	if atomCount > 2 and rOuterGroups < 1:
		temp = ""
		for symbol in rInnerParenthGroup:

			# Strip parenthesis in this >2 atom group, treat as linear group
			if symbol != '(':
				temp += symbol

			# Break at first atom so that expandGroup can finish the job instead of symbol by symbol logic inside data parser
			if symbol in ATOMS:
				break

		# Create new groups as symbol and first inner atom, bond included, and RN as the outer group with temp as inner group
		rInnerParenthGroupV1 = atomSymbol + '(' + temp + ')'
		rInnerParenthGroupV2 = atomSymbol + '(' + RN + ')' + temp
		# These parameters utilize the linearity of the SMILEScode to form a functional group from a ring,
		# using the two different fork paths that a ring junction creates
		# The form of the functionalGroup would be linear, as in RRR... with bonds and numbers.
		# it is created from the smiles representation of R(R...)R Where the inside of the parenthesis, the ...,
		# are removed. A functional groups connecting two branches of rings is thus created
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[-1], rInnerParenthGroupV1, [atomIndex, rInnerParenthindices[0]])
		if info is not False:
			for group in info:
				FGinfo.append(group)
		del(info)
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[0], rInnerParenthGroupV2, [atomIndex, RNindex, rInnerParenthGroup[0]])
		if info is not False:
			for group in info:
				FGinfo.append(group)
	elif atomCount == 1 and rOuterGroups < 1:

		rInnerParenthGroup = atomSymbol + rInnerParenthGroup
		# ##print("Evalutating rInnerParenthGroup", rInnerParenthGroup)
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[-1]+1, rInnerParenthGroup, [atomIndex, rInnerParenthindices[0]])
		if info is not False:
			for group in info:
				FGinfo.append(group)
	elif rOuterGroups == 2:
		# If two outer groups are attachted to a single atom, split the two groups by index and representaion
		rInnerPos = -1
		rInnerSmilesPos = rInnerParenthGroupPositions[0] - 2
		parenthesisCount = atomCount = 0
		firstGroup = True
		for symbol in rInnerParenthGroup:
			rInnerPos += 1
			rInnerSmilesPos += 1
			if symbol in ATOMS:
				atomCount += 1
			if symbol == '(':
				parenthesisCount += 1
			if symbol == ')':
				parenthesisCount -= 1
			if firstGroup is True and parenthesisCount == 0:
				firstEndParenthesis = rInnerPos + 1
				firstEndAtom = atomCount
				break
		firstParenthesisGroup = rInnerParenthGroup[0:firstEndParenthesis]
		firstParenthesisGroupindices = rInnerParenthindices[0:firstEndAtom]
		secondParenthesisGroup = '(' + rInnerParenthGroup[firstEndParenthesis:len(rInnerParenthGroup)]
		secondParenthesisGroupindices = rInnerParenthindices[firstEndAtom:len(rInnerParenthindices)]
		##print(firstParenthesisGroup, " with indices ", firstParenthesisGroupindices)
		##print(secondParenthesisGroup, " with indices ", secondParenthesisGroupindices)
		if firstParenthesisGroup != secondParenthesisGroup:
			rInnerParenthGroupV1 = atomSymbol + '(' + firstParenthesisGroup[1] + ')' + secondParenthesisGroup[1]
			info = whichGroup(atomPosition, rInnerSmilesPos + 1, rInnerParenthGroupV1, [atomIndex, firstParenthesisGroupindices[0], secondParenthesisGroupindices[0]])
			if info is not False:
				for group in info:
					FGinfo.append(group)
			rInnerParenthGroupV2 = atomSymbol + '(' + secondParenthesisGroup[1] + ')' + firstParenthesisGroup[1]
			info = whichGroup(atomPosition, rInnerSmilesPos + 1, rInnerParenthGroupV2, [atomIndex, secondParenthesisGroup[0], firstParenthesisGroup[0]])
			if info is not False:
				for group in info:
					FGinfo.append(group)
	elif rOuterGroups > 2:
		##print("There was an error")
		return 0

	# Combine RN info and determine whichGroup. If RN is empty, do not run this portion of code
	# ##print("RN = ", RN)
	if RN != "":
		RNgroup = atomSymbol + RNbond + RN
		#print("RNgroup is ", RNgroup)
		RNindices = [atomIndex, RNindex]
		info = whichGroup(atomPosition, RNPos, RNgroup, RNindices)
		if info is not False:
			for group in info:
				FGinfo.append(group)
		#print("found, ", info)


	# !!!!!!!!!!!!!  LN LOOP  !!!!!!!!!!!!!!!
	lOuterNumGroup = ""
	lBlockCounter = 0
	lOuterindices = []

	# Loop continues until an atom within the same grouping layer as atomSymbol is found to the left
	while LN not in ATOMS or lBlockCounter > 0:
		if LN in ATOMS:
			LNindex -= 1
		if LN == ')':
			lBlockCounter += 1
		if LN == '(' and lBlockCounter != 0:
			lBlockCounter -= 1
		if LN in BONDS and lBlockCounter == 0:
			LNbond = LN
		if LN in NUMBERS and lBlockCounter == 0:
			lOuterNumGroup = numbersHandler(LNPos)
		if LN in ATOMS and lBlockCounter != 0:
			lOuterindices.append(LNindex)
		LNPos -= 1
		LN = SMILEScode[LNPos]
	else:
		LNindex -= 1

	# Combine LN info and determine whichGroup
	if lOuterNumGroup:
		LNindices = [LNindex, lOuterNumGroup[1], atomIndex]
		LNgroup = LN + '(' +  lOuterNumGroup[0] + ')' + LNbond + atomSymbol
		#print("LNgroup with number is ", LNgroup)
	else:
		LNindices = [LNindex, atomIndex]
		LNgroup = LN + LNbond + atomSymbol

	# print("LN group is ", LNgroup, " with indices ", LNindices)
	info = whichGroup(LNPos, atomPosition, LNgroup, LNindices)
	# print("INFO = ",info)
	if info is not False:
		for group in info:
			FGinfo.append(group)

	# ##print("atomHandler found the group info to be, ", FGinfo)
	# If entries are in FGinfo, return it. Otherwise return False
	if FGinfo:
		return FGinfo
	else:
		return False


def chargeHandler(chargeSymbol, chargePostion, atomIndex):

	# Charges are a special case within SMILEScode. Each individual charge is enclosed
	# in brackets with the charged atom, in the form [C-] or [N+]. Charged hygroden
	# atoms are removed from the SMILEScode by the formatSMILEScode() function, since
	# no functioanl group needs to contain an expicitly charged hydrogen
	chargeGroup = SMILEScode[chargePostion-2:chargePostion+2]
	info = whichGroup(chargePostion-2, chargePostion+2, chargeGroup, [atomIndex])
	return info


def whichGroup(startPosition, endPosition, group, atomindices):


	# ##print("finding the group for ", group)
	# Holds templates which contain group, sent to expandGroup to attempt a fullMatch
	portionMatches = []

	# GITposition is cleared each time a whichGroup call is made to account for new group and template positional relation.
	GITposition.clear()

	# Only one fullMatch is possible
	fullMatch = False

	# Loop through groups
	for line in open('FGlist.txt', 'r'):

		# Variables
		lineInfo = re.compile(r'\S+').findall(line)
		FGtemplate = lineInfo[0].replace('[R]', 'R')
		FGname = lineInfo[1]
		difference = len(FGtemplate) - len(group)

		match = checkGroup(group, FGtemplate, False)
		# if FGname == "SecondaryAmine":
		# 	##print("match between ", group, " and ", FGtemplate, " is, ", match)
		# If group is larger than the template, or no match was found, skip
		if difference < 0 or match is False:
			continue

		# Equivalent lengths with a True match means an identical match
		elif difference == 0 and match is True:
			fullMatch = True
			fullDict = createFGDict(atomindices, FGtemplate)
			fullInfo = [FGname, FGtemplate, atomindices, fullDict]

		# Unequivalent lengths with a True match means a portion match
		elif difference > 0 and match is True:
			portionMatches.append([FGtemplate, FGname])

	# If a group was found without any portions, return group match information
	if fullMatch and not portionMatches:
		#print("fullMatch True, portionMatches False")
		return [fullInfo]

	# If a fullMatch and a portion was found, determine if on of the portions complete a group. Take expandGroup if it exists, otherwise take fullMatch
	elif fullMatch and portionMatches:
		expandedGroup = expandGroup(startPosition, endPosition, group, atomindices, portionMatches)
		if expandedGroup is not False:
			#print("Fullmatch with Expanded Group")
			return expandedGroup
		else:
			#print("Fullmatch with NO Expand Group")
			return [fullInfo]

	# If there are portions, determine if they complete a group
	elif fullMatch is False and portionMatches:
		expandedGroup = expandGroup(startPosition, endPosition, group, atomindices, portionMatches)
		if expandedGroup is not False:
			#print("No full match but Expand Group")
			return expandedGroup
		else:
			#print("No full match or expandGroup")
			return False

	# If there is no matches at all, return False
	elif fullMatch is False and not portionMatches:
		#print("No full match or portionMatches")
		return False


# portionMatches is a list of larger templates which contain the subtemplate, group, within it
# startPosition = string index of first symbol of a group within SMILEScode
# endPosition = string index of final symbol of a group within SMILEScode
# atomindices is a list of the indices which the group points to within the SMILEScode
def expandGroup(startPosition, endPosition, group, atomindices, portionMatches, recursive=False):

	# print("Expand group called with ", portionMatches)
	matches = []
	portionCounter = -1
	for portionGroup in portionMatches:
		portionCounter += 1
		# Initialize variables
		template = portionGroup[0]
		FGname = portionGroup[1]
		numTemplateAtoms = len(ATOMSREGEX.findall(template))
		difference = len(template) - len(group)
		posInTemplate = GITposition[portionCounter]
		requiredGroup = ""
		leftRequired = template[0:posInTemplate]
		rightRequired = template[posInTemplate+len(group):len(template)]
		numLeftRequired = len(leftRequired)
		numRightRequired = len(rightRequired)
		nonGroupTemplate = leftRequired + rightRequired
		numParenthesis = len(re.compile(r'\(').findall(nonGroupTemplate))
		# print("\nINITAL CONDITIONS ")
		# print("startPosition = ", startPosition)
		# print("endPosition = ", endPosition)
		# print("template = ", template)
		# print("group = ", group)
		# print("nonGroupTemplate = ", nonGroupTemplate)
		# print("FGname = ", FGname)
		# print("numTemplateAtoms = ", numTemplateAtoms)
		# print("posInTemplate = ", posInTemplate)
		# print("recursive = ", recursive)
		finalAtomindices = []
		if len(GITposition) != len(portionMatches) and len(portionMatches) != 1:
				print("GITposition error, the number of list entries is not equal to the number of portionMatches ")
				sys.exit()

		# If checkGroup finds a match, then left/rightExpand are false because the loop creation must be false
		# Likewise, if checkGroup does not find a match, then expansion logic is necessary
		leftExpansionFail = rightExpansionFail = False # Failure must be proven true in loop
		leftExpand = not checkGroup(group, template[0:len(group)], True)
		rightExpand = not checkGroup(group, template[posInTemplate:len(template)], True)
		if leftExpand is False:
			# this checks for a recursive call on rightExpand
			# Expand group can be rightward expanded in multiple directions
			# Mulitple expansion calls with various rightward variables passed
			# into expandGroup accounts for this. At this point in the expansion logic,
			# the leftward expansion must have been true already. This statement
			# checks if this is true by checking it the left expanded portion of the templates
			# has already been expanded into and is now held in the variable, group
			requiredGroup = group
			for index in atomindices:
				finalAtomindices.append(index)
		numRequiredGroupAtoms = len(ATOMSREGEX.findall(requiredGroup))
		if recursive is True and numRequiredGroupAtoms != len(atomindices):
			#print("atomindices IS popped")
			finalAtomindices.pop(-1)
			#print(atomindices)
		# print("leftExpand = ", leftExpand)
		# print("rightExpand = ", rightExpand)
		# print("requiredGroup = ", requiredGroup)
		# print("finalAtomindices = ", finalAtomindices)
		# print("atomindices = ", atomindices)
		# print("\n")
		# LN expansion

		# Varibles
		LNPos = startPosition
		LNIndex = atomindices[0]
		LNtempPos = posInTemplate

		# LN Loop

		# Loop can only be ceased with leftExpand True if all left symbols of the template are validaed in the SMILEScode
		# IF not, the group cannot exist and the right loop is not ran
		while leftExpand is True:
			LNPos -= 1
			LNtempPos -= 1
			#print("leftExpansionFail = ", leftExpansionFail)
			if len(requiredGroup) == numLeftRequired or leftExpansionFail is True:
				break
			if LNPos < 0:
				leftExpansionFail = True
				break
			LN = SMILEScode[LNPos]
			LNtemp = template[LNtempPos]
			# print("LNtemp = ", LNtemp)
			# print("LN = ", LN)
			# print("requiredGroup = ", requiredGroup)
			# six cases: LNtempplate can be in ATOMS, RGROUP, PARENTHESIS, CHARGES, BONDS, BRACKETS
			# If LNtemplate is in ATOMS, RGROUP, BONDS, BRACKET, and a PARENTHESIS or NUMBER is found, simply continue until the actual one is found
			if LNtemp == 'R':
				lBlockCounter = 0
				while LN not in ATOMS or lBlockCounter > 0:
					LNPos -= 1
					if LNPos < 0:
						leftExpansionFail = True
						break
					LN = SMILEScode[LNPos]
					if LN == ')':
						lBlockCounter += 1
					if LN == '(':
						lBlockCounter -= 1
					if lBlockCounter != 0 and LN in ATOMS:
						LNIndex -= 1
					if (lBlockCounter == 0 and LN not in ATOMS) or LNPos < 0:
						leftExpansionFail = True
						break
				else:
					LNIndex -= 1
					finalAtomindices.append(LNIndex)
					requiredGroup += LN
			elif LNtemp in BRACKETS:
				# Many molecules do not have any charges, so check if there are any  in the smilescode before executing the loop
				# If there are none, the template will never match. Break the loop set leftExpansionFail to False
				numberOfBrackets = len(re.compile(r'\[').findall(SMILEScode))
				if numberOfBrackets > 0:
					lBlockCounter = 0
					while LN not in BRACKETS or lBlockCounter > 0:
						LNPos -= 1
						if LNPos < 0:
							leftExpansionFail = True
							break
						LN = SMILEScode[LNPos]
						if LN == ')':
							lBlockCounter += 1
						if LN == '(':
							lBlockCounter -= 1
						if lBlockCounter != 0 and LN in ATOMS:
							LNIndex -= 1
						if (lBlockCounter == 0 and LN not in BRACKETS) or LNPos < 0:
							leftExpansionFail = True
							break
					else:
						requiredGroup += LN
				else:
					leftExpansionFail = True
					break
			elif LNtemp in ATOMS:
				lBlockCounter = 0
				while LN in NUMBERS or lBlockCounter > 0:
					LNPos -= 1
					if LNPos < 0:
						leftExpansionFail = True
						break
					LN = SMILEScode[LNPos]
					if LN == ')':
						lBlockCounter += 1
					if LN == '(':
						lBlockCounter -= 1
					if lBlockCounter != 0 and LN in ATOMS:
						LNIndex -= 1
					if LNPos < 0:
						leftExpansionFail = True
						break
				else:
					if LN == LNtemp:
						LNIndex -= 1
						finalAtomindices.append(LNIndex)
						requiredGroup += LN
					else:
						leftExpansionFail = True
						break
			elif LNtemp == ')':
				if LN in NUMBERS:
					numGroupinfo = numbersHandler(LNPos)
					numGroup = '(' + numGroupinfo[0] + ')'
					if numGroup == template[LNtempPos-len(numGroup)+1:LNtempPos+1]:
						requiredGroup += numGroup
					else:
						leftExpansionFail = True
						break
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
						LNPos -= 1
						if LNPos < 0:
							leftExpansionFail = True
							break
						LN = SMILEScode[LNPos]
						if LNPos in ATOMS:
							LNindex -= 1
						if LN == ')':
							lBlockCounter += 1
						if LN == '(':
							lBlockCounter -= 1
					else:
						outerGroup = SMILEScode[LNPos:LNPos+len(requiredGroup)-1] + ')'
						if outerGroup != requiredOuterGroup:
							leftExpansionFail = True
							break
						else:
							requiredGroup += outerGroup
				else:
					leftExpansionFail = True
					break
			elif LNtemp in BONDS:
				#print("LNtemp = ", LNtemp)
				#print("LN = ", LN)
				if LN in NUMBERS and (LNPos-1) >= 0:
					if SMILEScode[LNPos-1] == LNtemp:
						requiredGroup += SMILEScode[LNPos-1]
						LNPos -= 1
					else:
						leftExpansionFail = True
						break
				if LN == LNtemp:
					requiredGroup += LN
				else:
					leftExpansionFail = True
					break
			elif LNtemp in CHARGES:
				if LN not in CHARGES:
					leftExpansionFail = True
					break
				else:
					requiredGroup += LN
			elif LNtemp == '(':
				if LN != '(':
					leftExpansionFail = True
					break
				else:
					requiredGroup += LN
		if leftExpansionFail is True: # If loop failed, then match is impossible within SMILEScode
			#print("leftExpansion failed, continuing")
			continue
		# If there was a leftExpansion and the loop did not fail, then attatch the group and indices to finalAtomindices
		# If there was no leftExpansion, requiredGroup is already equal to group, see above the loop
		elif leftExpansionFail is False and leftExpand is True:
			finalAtomindices = finalAtomindices[::-1]
			requiredGroup = requiredGroup[::-1]
			requiredGroup += group
			for index in atomindices:
				finalAtomindices.append(index)

		# RN expansion

		# Variables

		RNPos = endPosition # RNPosition in SMILEScode
		RNIndex = atomindices[-1] # RNindex in SMILEScode
		RNtempPos = len(requiredGroup) - 1 # RNtemplateposition, intialized at final position of group within template
		if leftExpand is False:
			RNtempPos = len(group) - 1
		while rightExpand is True:
			RNPos += 1
			RNtempPos += 1
			if len(requiredGroup) == len(template): # If requiredGroup reaches the length of the template, then a match has been found
				# print(requiredGroup, " == ", template)
				break
			if RNPos > SMILEScodelength - 1 or rightExpansionFail is True:
				rightExpansionFail = True
				break
			RN = SMILEScode[RNPos]
			RNtemp = template[RNtempPos]
			# print("RNtemp is = ", RNtemp)
			# print("RN = ", RN)
			# print("requiredGroup = ", requiredGroup)
			# print("RNtempPos = ", RNtempPos)
			if RN == '(':
				rBlockCounter = 1
				correlatingParenthPos = RNPos
				correlatingParenth = SMILEScode[RNPos]
				outerIndex = RNIndex
				innerInidicies = [outerIndex]
				while rBlockCounter > 0:
					correlatingParenthPos += 1
					correlatingParenth = SMILEScode[correlatingParenthPos]
					if correlatingParenth in ATOMS:
						outerIndex += 1
						innerInidicies.append(outerIndex)
					if correlatingParenth == '(':
						rBlockCounter += 1
					if correlatingParenth == ')':
						rBlockCounter -= 1
				innerSymbol = SMILEScode[RNPos+1]
				innerIndex = RNIndex + 1
				outerSymbol = SMILEScode[correlatingParenthPos+1]
				outerIndex += 1
			if RNtemp == 'R':
				if RN == '(':
					# Only an Rgroup can exist after the parenthesis
					if innerSymbol in ATOMS:
						InnerExpansion = expandGroup(LNPos, RNPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if InnerExpansion is not False and recursive is False:
							matches.append(InnerExpansion)
							rightExpansionFail = True
						elif InnerExpansion is not False and recursive is True:
							return InnerExpansion
					elif outerSymbol in ATOMS:
						OuterExpansion = expandGroup(LNPos, correlatingParenthPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if OuterExpansion is not False and recursive is False:
							matches.append(OuterExpansion)
							rightExpansionFail = True
						elif OuterExpansion is not False and recursive is True:
							return OuterExpansion
				if RN in NUMBERS:
					numGroup = numbersHandler(RNPos)
					if numGroup[3] == "Linear":
						requiredGroup += numGroup[0]
						finalAtomindices.append(numGroup[1])
					elif (RNPos+1) < SMILEScodelength:
						if SMILEScode[RNPos+1] in ATOMS:
							RNIndex += 1
							requiredGroup += SMILEScode[RNPos+1]
							RNPos += 1
							finalAtomindices.append(RNIndex)
						else:
							rightExpansionFail = True
							break
					else:
						rightExpansionFail = True
						break
				elif RN in ATOMS:
					RNIndex += 1
					requiredGroup += RN
					finalAtomindices.append(RNIndex)
				else:
					rightExpansionFail = True
					break
			elif RNtemp in ATOMS:
				if RN == '(':
					# The specific atom can exist only after the parenthesis
					if innerSymbol == RNtemp:
						InnerExpansion = expandGroup(LNPos, RNPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if InnerExpansion is not False and recursive is False:
							matches.append(InnerExpansion)
							rightExpansionFail = True
						elif InnerExpansion is not False and recursive is True:
							return InnerExpansion
					elif outerSymbol == RNtemp:
						OuterExpansion = expandGroup(LNPos, correlatingParenthPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if OuterExpansion is not False and recursive is False:
							matches.append(OuterExpansion)
							rightExpansionFail = True
						elif OuterExpansion is not False and recursive is True:
							return OuterExpansion
				if RN in NUMBERS:
					numGroup = numbersHandler(RNPos)
					if numGroup[3] == "Linear" and numGroup[0] == RNtemp:
						requiredGroup += numGroup[0]
						finalAtomindices.append(numGroup[1])
					elif (RNPos+1) < SMILEScodelength:
						if SMILEScode[RNPos+1] == RNtemp:
							RNIndex += 1
							requiredGroup += SMILEScode[RNPos+1]
							RNPos += 1
							finalAtomindices.append(RNIndex)
						else:
							rightExpansionFail = True
							break
					else:
						rightExpansionFail = True
						break
				elif RN == RNtemp:
					RNIndex += 1
					requiredGroup += RN
					finalAtomindices.append(RNIndex)
				else:
					rightExpansionFail = True
					break
			elif RNtemp in BRACKETS:
				if RN == '(':
					# Only an Rgroup can exist after the parenthesis
					if innerSymbol in BRACKETS:
						InnerExpansion = expandGroup(LNPos, RNPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if InnerExpansion is not False and recursive is False:
							matches.append(InnerExpansion)
							rightExpansionFail = True
						elif InnerExpansion is not False and recursive is True:
							return InnerExpansion
					elif outerSymbol in BRACKETS:
						OuterExpansion = expandGroup(LNPos, correlatingParenthPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if OuterExpansion is not False and recursive is False:
							matches.append(OuterExpansion)
							rightExpansionFail = True
						elif OuterExpansion is not False and recursive is True:
							return OuterExpansion
				if RN in NUMBERS:
					numGroup = numbersHandler(RNPos)
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
						else:
							rightExpansionFail = True
							break
					elif (RNPos+1) < SMILEScodelength:
						if SMILEScode[RNPos+1] in BRACKETS:
							requiredGroup += RN
							RNPos += 1
						else:
							rightExpansionFail = True
							break
					else:
						rightExpansionFail = True
						break
				elif RN in BRACKETS:
					requiredGroup += RN
				else:
					rightExpansionFail = True
					break
			elif RNtemp in BONDS:
				if RN == '(':
					# Only an Rgroup can exist after the parenthesis
					if innerSymbol in BONDS:
						InnerExpansion = expandGroup(LNPos, RNPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if InnerExpansion is not False and recursive is False:
							matches.append(InnerExpansion)
							rightExpansionFail = True
						elif InnerExpansion is not False and recursive is True:
							return InnerExpansion
					elif outerSymbol in BONDS:
						OuterExpansion = expandGroup(LNPos, correlatingParenthPos, requiredGroup, finalAtomindices, [portionGroup], True)
						if OuterExpansion is not False and recursive is False:
							matches.append(OuterExpansion)
							rightExpansionFail = True
						elif OuterExpansion is not False and recursive is True:
							return OuterExpansion
				if RN in NUMBERS and (RNPos + 1) < SMILEScodelength:
					if SMILEScode[RNPos+1] == RNtemp:
						requiredGroup += SMILEScode[RNPos+1]
						RNPos += 1
					else:
						rightExpansionFail = True
						break
				elif RN == RNtemp:
					requiredGroup += RN
				else:
					rightExpansionFail = True
					break
			elif RNtemp in CHARGES:
				if RN in CHARGES:
					requiredGroup += RN
				else:
					rightExpansionFail = True
					break
			elif RNtemp == '(':
				parenthesiedTemp = '('
				while RNtemp != ')':
					RNtempPos += 1
					RNtemp = template[RNtempPos]
					parenthesiedTemp += RNtemp
				# print("parenthesiedTemp = ", parenthesiedTemp)
				if RN in NUMBERS:
					numGroup = numbersHandler(RNPos)
					# If group is outer, add it. If it is linear, incompatable and group does not exist where number is represented at
					if numGroup[3] == "Outer":
						outerParenthesied = '(' + numGroup[0] + ')'
					else:
						rightExpansionFail = True
						break
					validNumGroup = checkGroup(outerParenthesied, parenthesiedTemp, True)
					# print("validNumGroup = ", validNumGroup)
					if validNumGroup is False:
						rightExpansionFail = True
						break
					else:
						requiredGroup += outerParenthesied
						finalAtomindices.append(numGroup[1])
						# RNtempPos += (len(outerParenthesied) - 2s)
				elif RN == '(':
					innerParenthesied = SMILEScode[RNPos:RNPos+len(parenthesiedTemp)-1] + ')'
					outerParenthesied = '(' + SMILEScode[correlatingParenthPos+1:correlatingParenthPos+1+len(parenthesiedTemp)-2] + ')'
					innerExpand = checkGroup(innerParenthesied, parenthesiedTemp, True)
					outerExpand = checkGroup(outerParenthesied, parenthesiedTemp, True)
					# print("innerParenthesied = ", innerParenthesied)
					# print("outerParenthesied = ", outerParenthesied)
					if innerExpand is True:
						indices = []
						for index in finalAtomindices:
							indices.append(index)
						numAtoms = len(ATOMSREGEX.findall(innerParenthesied))
						for number in range(0,numAtoms):
							indices.append(innerIndex+number)
						indices.append(outerIndex-1) # Ensure that OuterExpansion begins on the correct index, see pop() call on line 857
						# print("Outer expanding with ", requiredGroup+innerParenthesied, " and the index's", indices)
						OuterExpansion = expandGroup(startPosition, correlatingParenthPos, requiredGroup+innerParenthesied, indices, [portionGroup], True)
						if OuterExpansion is not False and recursive is False:
							matches.append(OuterExpansion)
							rightExpansionFail = True
						elif OuterExpansion is not False and recursive is True:
							return OuterExpansion
					elif outerExpand is True:
						indices = []
						for index in finalAtomindices:
							indices.append(index)
						numAtoms = len(ATOMSREGEX.findall(outerParenthesied))
						for number in range(0,numAtoms):
							indices.append(outerIndex+number)
						indices.append(innerIndex-1) # Ensure that InnerExpansion begins on the correct index, see pop() call on line 857
						# print("Inner expanding with ", requiredGroup+ outerParenthesied, " and the indices ", indices)
						InnerExpansion = expandGroup(startPosition, RNPos, requiredGroup+outerParenthesied, indices, [portionGroup], True)
						if InnerExpansion is not False and recursive is False:
							matches.append(InnerExpansion)
							rightExpansionFail = True
						elif InnerExpansion is not False and recursive is True:
							return InnerExpansion
					if outerParenthesied != parenthesiedTemp and innerParenthesied != parenthesiedTemp:
						rightExpansionFail = True
						break
				else:
					rightExpansionFail = True
					break
			elif RNtemp == ')':
				if RN in NUMBERS and (RNPos+1) < SMILEScodelength:
					if SMILEScode[RNPos+1] == ')':
						requiredGroup += SMILEScode[RNPos+1]
						RNPos+=1
					else:
						rightExpansionFail = True
						break
				elif RN == ')':
					requiredGroup += RN
				else:
					rightExpansionFail = True
					break

		if rightExpansionFail is True:
			#print("rightExpansionFail = ", rightExpansionFail)
			continue
		elif rightExpansionFail is False:
			# print("rightExpansionFail = ", rightExpansionFail)
			finalCheck = checkGroup(requiredGroup, template, True)
			# print("finalCheck = ", finalCheck)
			# print(len(finalAtomindices), " == ", numTemplateAtoms)
			# A final check for if the group is identical to requiredGroup and the number of atom indices is equal to the number of atoms in the template
			if finalCheck is True and len(finalAtomindices) == numTemplateAtoms:
				finalDict = createFGDict(finalAtomindices, template)
				FGinfo = [FGname, template, finalAtomindices, finalDict]
				# print(FGinfo)
				# print(len(portionMatches))
				if len(portionMatches) == 1 and recursive is True: # If only one portionMatch, then it must be a recursive call because expandGroup always has more than one portion
					# print("Returning SINGLE FGinfo for RECURSIVE CALL\nrecursive = ", recursive)
					# print("FGInfo = ", FGinfo)
					# print("RECURISVE END!")
					return FGinfo
				else:
					matches.append(FGinfo)
			else:
				print("Fatal error, both expansion logics were true but the final check was false")
				print("finalCheck = ", finalCheck)
				print("finalAtomindices = ", finalAtomindices)
				print("numTemplateAtoms = ", numTemplateAtoms)
				print("group = ", group)
				print("requiredGroup = ", requiredGroup)
				print("template = ", template)
				sys.exit()

	# print("Returning matches = ", matches)
	if matches:
			# print(group)
			# for group in matches:
		return matches
	else:
		return False


def numbersHandler(numPosition):

	for RINGGROUP in RINGPOSITIONS:

		# If the number OPENS the group, return the CLOSING info
		if RINGGROUP[0][0] == numPosition:
			return(RINGGROUP[1][1], RINGGROUP[1][2], RINGGROUP[1][3], RINGGROUP[0][4])

		# If the number CLOSES the group, return the OPENING info
		elif RINGGROUP[1][0] == numPosition:
			return(RINGGROUP[0][1], RINGGROUP[0][2], RINGGROUP[0][3], RINGGROUP[1][4])


def createFGDict(FGindices, FGtemplate):
	FGdict = []
	indexCounter = -1
	#print(FGindices)
	numAtoms = len(ATOMSREGEX.findall(FGtemplate))
	#print(numAtoms)
	for symbol in FGtemplate:
		if symbol in ATOMS or symbol == 'R':
			indexCounter += 1
			#print(indexCounter)
			FGdict.append([symbol, FGindices[indexCounter]])
	#print(FGdict)
	return FGdict


def checkGroup(group, template, full=True):

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
			if groupSymbol not in RGROUP and symbol == 'R': # If symbol is an atom where it should be
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
				if symbol not in RGROUP and templateSymbol == 'R':
					break
				if symbolCounter == len(group) - 1:
					GITposition.append(shift)
					return True
		return False


def formatSMILEScode(SMILEScode):
	SMILEScodePos = -1
	reFormatted = ""
	for symbol in SMILEScode:
		SMILEScodePos += 1
		if symbol == '[':
			startBracketPos = SMILEScodePos
		if symbol == 'H':
			cutPos = SMILEScodePos
			while SMILEScode[cutPos] != '+':
				cutPos += 1
			reFormatted = SMILEScode[0:startBracketPos] + SMILEScode[startBracketPos+1] + SMILEScode[cutPos+2:len(SMILEScode)]
			reFormatted = formatSMILEScode(reFormatted)
			break
	if reFormatted == "":
		return SMILEScode
	else:
		return reFormatted


def initializeCYCLICINDICES():

	evaluatedNumbers = []

	SMILEScodePos = atomIndex = -1
	for symbol in SMILEScode:
		SMILEScodePos+=1

		if symbol in ATOMS:
			atomIndex += 1

		if symbol in NUMBERS and SMILEScodePos not in evaluatedNumbers:
			# print("Found ", symbol, " at position, ", SMILEScodePos)
			# print("\n")

			outerIndices = []
			rBlockCounter = 0 # Counter to ensure same scope of indicies is appended to CYCLICINDICES
			evaluatedNumbers.append(SMILEScodePos)
			CYCLICINDICES.append(atomIndex)
			# print(SMILEScode[SMILEScodePos])
			RNPos = SMILEScodePos + 1
			RNindex = atomIndex
			# print(SMILEScode)
			RN = SMILEScode[RNPos]

			while RN != symbol:
				# print("RN = ", RN)
				# print("rBlockCounter = ", rBlockCounter)
				if RN == '(':
					rBlockCounter += 1
				if RN == ')':
					rBlockCounter -= 1
				if rBlockCounter == 0:
					outerIndices.clear() # Refresh every scope reset
				if RN in ATOMS:
					RNindex+=1
					if rBlockCounter == 0:
						CYCLICINDICES.append(RNindex)
					elif rBlockCounter != 0:
						outerIndices.append(RNindex)
				RNPos+=1
				RN = SMILEScode[RNPos]
			else:
				# print("Correlating number", RN, " found at ", RNPos)
				evaluatedNumbers.append(RNPos)
				if rBlockCounter != 0:
					for index in outerIndices:
						CYCLICINDICES.append(index)

	# Distinguish between aromatic and non aromatic indicies
	SMILEScodePos = atomIndex = -1
	for symbol in SMILEScode:
		SMILEScodePos += 1
		if symbol in ATOMS:
			atomIndex+=1
			if symbol.islower() and atomIndex not in AROMATICINDICIES:
				AROMATICINDICIES.append(atomIndex)

def determineCyclicGroups(FGdataFinal):
	groupCounter = -1
	for group in FGdataFinal:
		groupCounter += 1
		groupIndices = group[2]

		# Check for Aromatic first, since aromatic group implies cyclic group
		isAromatic = all(index in AROMATICINDICIES for index in groupIndices)
		if isAromatic:
			group[0] = "Aromatic" + group[0]
			continue

		# Check for case of all indices being cyclic, automatically a cyclic group
		isCyclic = all(index in CYCLICINDICES for index in groupIndices)
		if isCyclic:
			group[0] = "Cyclic" + group[0]
			continue

		# Check for case when at least 2 consecutive indices are cyclic, then group is cyclic
		sampleIndices = []
		indexCounter = -1
		while isCyclic is False:
			indexCounter += 1
			# Checking one ahead, and 0 index of lists. So break at length - 2
			if indexCounter == len(groupIndices) - 2:
				break
			sampleIndices.append(groupIndices[indexCounter])
			sampleIndices.append(groupIndices[indexCounter+1])
			isCyclic = all(index in CYCLICINDICES for index in sampleIndices)
			sampleIndices.clear() # Refresh for a new set each increment of indexCounter
		if isCyclic:
			group[0] = "Cyclic" + group[0]
			continue
