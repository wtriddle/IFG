# Written by William Riddle from 9/23/19 -

import rdkit.Chem
import re
import collections


# The data parsing functions written here slice the SMILEScode in various ways to obtain a string of a potential functionalGroup, represented in the FGlist.txt
# To optimze the program, the he data parsing functions can find the same group multiple times, such that a single functionalgroup may appear multiple times in FGdata
# This is intended so that the SMILEScode can be looked at from various viewpoints, using its linearity

# VARIBLE DESCRIPTIONS

# RN stands for rightNeighbor, and evalutes symbols to the right of the symbol of interest in the SMILES code
# LN stands for leftNeighbor, and evalutes symbols to the left of the symbol of interest in the SMILES code

# FGinfo = (FGname, FGtemplate, FGindicies, FGdict) format of extracted functional group data
# FGname = Name of Functional group
# FGtemplate = SMILEScode equivalent of functional group i.e. C(=O) is template for ketone
# FGindicies = atomIndicies involved in group
# FGdict = dictionary corellating indicies to template positions. In the form for n atoms in a group, ['templateSymbol', index]. i.e. ['R', 3] or ['n', 6]

# RINGPOSITIONS = global variable which holds the opening and closing atom data for a given ring. Supersedes linearity of SMILEScode to implement the proper dynamic ring structure
# to model the chemical strucutre in space. String of characters -> Dynamic strucutre of atoms.
# For (2n - n) numbers in the SMILEScode, RINGPOSITIONS has n cells in the form:
# ([openingPosition, openingIndex, openingSymbol], [closingPosition, closingIndex, closingSymbol])

# FGdata = global variable containing every functional group identified by all data parsing functions
# For n number of FG's identified, FGdata has the cell form:
# [FGnfo]

RINGPOSITIONS = []
FGdata = []

# Constants
ATOMS = ['C','H','O','N','c','h','n','o']
RGROUP = ['C','H','O','N','c','h','n','o']
CHARGES = ['+','-']
PARENTHESIS = ['(', ')']
BRACKETS = ['[', ']']
BONDS = ['=', '#']
NUMBERS = ['1','2','3','4','5','6','7','8','9']


# Head function to loop through SMILEScode and call proper data parsers on given characters. Final FGdata is formatted as well
def ifg(SMILES, mol):
	global SMILEScode = SMILES  # Definition of SMILEScode
	global SMILESlength = len(SMILES)
	atomIndex = smilesPos = -1
	for symbol in SMILEScode:
		smilesPos += 1

		if symbol in PARENTHESIS:
			continue

		elif symbol in CHARGES:
			chargeGroup = chargeHandler(symbol, smilesPos, atomIndex)
			if chargeGroup is not FALSE:
				for group in chargeGroup:
					FGdata.append(chargeGroup)
			del(chargeGroup)

		elif symbol in ATOMS and smilesPos != 0 and smilesPos != SMILESlength - 1:
			atomIndex += 1
			directGroup = atomHandler(symbol, smilesPos, atomIndex)
			if directGroup is not FALSE:
				for group in directGroup:
					FGdata.append(directGroup)
			del(directGroup)

		elif symbol in BONDS:
			bondGroup = bondHandler(symbol, smilesPos, atomIndex)
			if bondGroup is not FALSE:
				for group in bondGroup:
					FGdata.append(bondGroup)
			del(bondGroup)

	FGdataFinal = evaluateFGdata()
	return FGdataFinal


# Initializes the global RINGPOSITIONS according to its model
def initializeRINGPOSITIONS():

	# Initialize variables
	smilesPos = atomIndex =  -1
	evaluatedNumbers = []
	chargeGroup = ""

	# Main loop
	for symbol in SMILEScode:
		smilesPos += 1

		# Keep track of ATOM positional data for when a number is found
		if symbol in ATOMS:
			atomIndex += 1
			# Atom found before number is opening atom of the ring.
			openingAtom = symbol
			openingAtomPos = smilesPos

		# When an unevalutaed number is found, evalute it
		if symbol in NUMBERS and smilesPos not in evaluatedNumbers:
			# smilesPos pointing at the position of the atom symbol in SMILEScode 1 space before the first occurance of the number

			# However, if it is not an atom in this scope, then it must be a bracket surrounding a charged atom.
			# The following if statement below is a special case where charged atom is opening the ring.
			# A closing bracket must precede the position of the number.
			if SMILEScode[smilesPos-1] == ']' and smilesPos >= openingAtomPos + 2:
				specialCounter = smilesPos - 1 # Initialized at posion of bracket
				specialSymbol = SMILEScode[specialCounter]
				while specialSymbol != '[':
					chargeGroup.insert(0, specialSymbol)
					specialCounter -= 1
					specialSymbol = SMILEScode[smilesPos]
				chargeGroup.insert(0, '[')

				# Reassign openingAtom and openingAtomPos so that the data can be \
				# appended to RINGPOSITIONS with a single statement
				openingAtom = chargeGroup
				openingAtomPos = specialCounter

			# subLoop counters
			correlatingPos = smilesPos - 1
			correlatingIndex = atomIndex

			# subLoop to find correlating atom in ring
			for subSymbol in SMILEScode[smilesPos:len(SMILEScode)]:

				correlatingPos += 1

				# Keep track of ATOM positional data for when the correlating atom is found
				if subSymbol in ATOMS:
					correlatingAtom = subSymbol
					correlatingIndex += 1
					correlatingAtomPos = correlatingPos

				# If the same number which began the sub loop is found, then append both opening and closing info
				if subSymbol == symbol:

					# If the ring opens with a charge, it must close with a charge. Re-assign positional variables
					if chargeGroup != "":
						correlatingChargeGroup = ""
						specialCounter = correlatingPos -1
						specialSymbol = SMILEScode[correlatingPos]
						while specialSymbol != '[':
							correlatingChargeGroup.insert(0, specialSymbol)
							specialCounter -= 1
							specialSymbol = SMILEScode[specialCounter]
						correlatingAtom = correlatingChargeGroup
						correlatingAtomPos = specialCounter

					RINGPOSITIONS.append([
						smilesPos,
						openingAtom,
						atomIndex,
						openingAtomPos],
						[correlatingPos,
						correlatingAtom,
						correlatingIndex,
						correlatingAtomPos])

					# Add numbers to evaluted list to prevent reevalutaion of correlatingAtom in main loop
					evaluatedNumbers.append(smilesPos)
					evaluatedNumbers.append(correlatingPos)
					break

				# Error caught if a correlatingPos runs out of the scope of the SMILEScode
				if correlatingPos == len(SMILEScode) - 1:
					print("There has been an error in the creation of the RINGPOSITIONS global used in the program\n")
					print("The smiles code ", SMILEScode, " had an error")
	return 0


# evaluateFGdata() evaluates the global FGdata variable and determines the valid groups within the strucutre, with proper occurances
# The formattedFGdata list takes from the FGdata data under certain contions. The following are descrbed
# Given a group and a compareGroup, both within FGdata:
#
#
# Note that the data parsers can label the same FG multiple times, (1) and (2) ensure this does not cause invalid data
# by checking if an identical group evaluted to be valid has already been validated or is comparing itself.
# A. group will be skipped if:
# (1) group is already within formattedFGdata
# (2) group is equal to compareGroup (i.e. FGdata[3] == FGdata[3])
#
# B and C are the group containment statements
# B. compareGroup will be skipped if:
# (1) compareGroup is fully contained within group (i.e Ketone contained in Ester)
#
# C. compareGroup will be added if:
# (1) group is fully contained within compareGroup, and compareGroup is not within formattedFGdata
#
# If no comparitive conclusion is made, then it is an individual group found by the data parsers.
# D. group will be added if:
# (1) compareGroup above was not added
# (2) group is not already within formattedFGdata
def evaluateFGdata():
	formattedFGdata = []
	groupCounter = -1
	for group in FGdata:

		groupCounter += 1

		# group variables
		groupTemplate = group[1]
		groupIndicies = group[2]
		groupDict = group[3]
		compareCounter = -1

		# Check if group is already in formattedFGdata (A.1) Skip if true
		inFormattedData = isGroupInList(formattedFGdata, group)
		if inFormattedData is TRUE:
			continue

		for compareGroup in FGdata:

			# compareGroup variables
			compareCounter += 1
			compareTemplate = compareGroup[1]
			compareIndicies = compareGroup[2]
			compareDict = compareGroup[3]

			# Check if loop is comparing the exact same group (A.2) Skip if true
			if compareCounter == groupCounter:
				continue

			# Determine number of Maingroup atoms
			numMainGroupAtoms = len(re.compile(r'R').findall(groupTemplate)) - len(groupIndicies)
			numMainCompareAtoms = len(re.compile(r'R').findall(compareTemplate)) - len(compareIndicies)

			# Determine Maingroup crossover
			mainGroupCrossover = 0
			for atom in groupDict:
				for compareAtom in compareDict:
					if atom[0] == compareAtom[0] and atom[0] != 'R' and compareAtom[0] != R: # If the symbols are symbol equivalent and both main group Atoms
						if atom[1] == compareAtom[1]: # If the atoms are the same index
							mainGroupCrossover += 1 # Then atoms must be of the exact same maingroup in the structure

			# Check if compareGroup is fully contained within group (B.1) If ture, skip this compareGroup.
			if mainGroupCrossover == numMainCompareAtoms and len(groupIndicies) > len(compareIndicies):
				continue

			# Check if group is fully contained within compareGroup (C.1) If so, and is a new group, add compareCounter
			if mainGroupCrossover == numMainGroupAtoms and len(groupIndicies) < len(compareIndicies):
				inFormattedData = isGroupInList(formattedFGdata, compareGroup)
				if inFormattedData is FALSE:
					formattedFGdata.append(compareGroup)
					break

		# If the group is not already in formattedFGdata, and it is not contained within another group (D.1, D.2) Add the group
		formattedFGdata.append(group)

	return formattedFGdata


# Called from evaluateFGdata, checks if a group is already contained within the list.
# Only EXACT matches, i.e fginfo === fginfo for all cells is true
def isGroupInList(formattedFGdata, group):

	# If there is an atom in formattedFGdata
	if formattedFGdata:

		for compareGroup in formattedFGdata:

			# Comparative Variables
			compareGroupTemplate = compareGroup[1]
			compareGroupIndicies = compareGroup[2]
			compareGroupDict = compareGroup[3]
			numMainCompareAtoms = len(re.compile(r'R').findall(compareGroupTemplate)) - len(compareGroupTemplate)

			# Group Variables
			groupTemplate = group[1]
			groupIndicies = group[2]
			groupDict = group[3]
			numMainGroupAtoms = len(re.compile(r'R').findall(group)) - len(groupTemplate)

				# Check for equivalent number of atoms
				if len(compareGroupIndicies) == len(groupIndicies):
					mainGroupCrossover = 0
					for atom in groupDict:
						for compareAtom in compareGroupDict:

							# Check if atoms are maingroup equivalent, and if they are in the same positional index.
							# The atoms in both functionalGroup scopes are "identical" if so, so add to mainGroupCrossover
							if atom[0] == compareAtom[0] and atom[0] != 'R':
								if atom[1] == compareAtom[1]:
									mainGroupCrossover += 1

					# If there is full maingroups crossover and the templates are the same
					if mainGroupCrossover == numMainCompareAtoms == numMainGroupAtoms and compareGroupTemplate == groupTemplate:
							# Then the exact group is already contained within formattedFGdata
							return TRUE
	return FALSE


# Creates groups stemming from BONDS symbol within SMILEScode
def bondHandler(bondSymbol, bondPostion, atomIndex):

	# Initialize positional variables to create group
	LNPos = bondPostion - 1
	RNPos = bondPostion + 1
	LNindex = RNindex = atomIndex
	LN = SMILES[LNPos]
	RN = SMILES[RNPos]


	# Special case which occurs commonly as R(=R)
	if LN == '(' and SMILES[LNPos-1] in ATOMS and RN in ATOMS and SMILES[RNPos+1] == ')':
		atomIndicies = [atomIndex-1, atomIndex]
		bondGroup = SMILES[LNPos-1] + "(=" + LN + ')'
		info = whichGroup(RNPos+1, LNPos-1, bondGroup, atomIndicies)
		if info is not FALSE:
			return info
		else:
			return FALSE


	# LN loop
	lBlockCounter = 0
	lOuterParenthGroup = ""
	lOuterNumGroup = ""


	# Find the first non blocking group atom connected to the bond
	# If LN is an atom, loop not run
	while LN not in ATOMS and lBlockCounter != 0:
		LNPos -= 1
		LN = SMILES[LNPos]
		if LN == ')':
			lBlockCounter += 1
		if LN == '(' and lBlockCounter != 0:
			lBlockCounter -= 1
		if LN in NUMBERS and lBlockCounter == 0:
			lOuterNumGroup = numbersHandler(LNPos)
		if lBlockCounter != 0:
			lOuterParenthGroup.insert(0, symbol)
		if LN in ATOMS:
			LNindex -= 1

	# RN loop. If RN is already an atom, loop not run
	# The atom next to the bond symbol is always the atom involved in the bond, no blocking groups occur
	while RN not in ATOMS:
		RNPos += 1
		RN = SMILES[RNPos]
		if RN in ATOMS:
			RNindex += 1

	# lOuterParenthGroup evaluation to see if the group should be included in bondGroup
	atomCount = len(re.compile(ATOMS).findall(lOuterParenthGroup))
	numCount = len(re.compile(NUMBERS).findall(lOuterParenthGroup))
	if atomCount > 2 or numCount != 0:
		lOuterParenthGroup = ""

	# Bond Group creation
	if lOuterNumGroup:
		atomIndicies = [LNindex, lOuterNumGroup[1], RNindex]
		bondGroup = LN + lOuterNumGroup[0] + lOuterParenthGroup + bondSymbol + RN
	else:
		atomIndicies = [LNindex, RNindex]
		bondGroup = LN + lOuterParenthGroup + bondSymbol + RN

	# Determine the group and return the info
	FGinfo = whichGroup(RNPos, LNPos, bondGroup, atomIndicies)
	return FGinfo



# Creates groups stemming from ATOMS symbol within SMILEScode
# Special RN case for first inner atom + outer atom group not written yet. Place in if needed.
def atomHandler(atomSymbol, atomPosition, atomIndex):

	# Initialize Variables
	LNPos = atomPosition - 1
	RNPos = atomPosition + 1
	RN = SMILES[RNPos]
	LNPos = SMILES[LNPos]
	RNindex = LNindex = atomIndex
	RNbond = LNbond = ""
	FGinfo = []

	# ----RN LOOP----
	rBlockCounter = 0
	rInnerParenthGroup = ""
	rInnerParenthIndicies = []
	rInnerParenthGroupPositions = []

	# Loop continues until an atom within the same grouping layer as atomSymbol is found to the right
	while RN not in ATOMS and rBlockCounter != 0:

		# IF RN in the same grouping layer as atomSymbol
		RNPos += 1
		RN = SMILES[RNPos]
		if RN in BONDS and rBlockCounter == 0:
			RNbond = RN
		if RN in ATOMS:
			RNindex += 1
		if RN == '(':
			rBlockCounter += 1
		if RN == ')':
			rBlockCounter -= 1

		# If RN on the edge of a parenthesis, no RN exists, cease the loop
		if rBlockCounter < 0:
			RN = ""
			break

		# If scope of RN loop is within a parenthesis, capture the information inside of it with respect to atomSymbol
		if rBlockCounter != 0:
			rInnerParenthGroup += RN
			rInnerParenthGroupPositions.append(RNPos)
		if rBlockCounter != 0 and RN in ATOMS:
			rInnerParenthIndicies.append(RNindex)
		if rBlockCounter != 0 and RN in NUMBERS:
			numGroup = numbersHandler(RNPos)
			rInnerParenthIndicies.append(numGroup[1])
			rInnerParenthGroup += numGroup[0]


	# rInnerParenthGroup evaluation
	atomCount = len(re.compile(ATOMS).findall(rInnerParenthGroup))

	if atomCount > 2:
		temp = ""
		for symbol in rInnerParenthGroup:

			# Strip parenthesis in this >2 atom group, treat as linear group
			if symbol != '(':
				temp = symbol

			# Break at first atom so that expandGroup can finish the job instead of symbol by symbol logic inside data parser
			if symbol in ATOMS:
				break

		# Create new group as symbol and first inner atom, bond included
		rInnerParenthGroup = atomSymbol + temp

		# These parameters utilize the linearity of the SMILEScode to form a functional group from a ring,
		# using the two different fork paths that a ring junction creates
		# endPosition in whichGroup, parameter 2, is the final position of rInnerParenthGroup.
		# Since the rInnerParenthGroupPositions tracks the position, the last element is passed in.
		# The added 1 is to compensate for the final ')', to begin at the other forked path
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[-1] + 1, rInnerParenthGroup, rInnerParenthIndicies[0])
		if info is not FALSE:
			FGinfo.append(info)

	elif atomCount == 1:

		# Add final ')' to group
		rInnerParenthGroup = atomSymbol + rInnerParenthGroup + ')'
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[0], rInnerParenthGroup, rInnerParenthIndicies[0])
		if info is not FALSE:
			FGinfo.append(info)

	# Combine RN info and determine whichGroup
	RNgroup = atomSymbol + RNbond + RN
	RNindicies = [atomIndex, RNindex]
	info = whichGroup(atomPosition, RNPos, RNgroup, RNindicies)
	if info is not FALSE:
		FGinfo.append(info)


	# ----LN LOOP----
	lOuterNumGroup = ""
	lBlockCounter = 0

	# Loop continues until an atom within the same grouping layer as atomSymbol is found to the left
	while LN not in ATOMS and lBlockCounter != 0:

		LNPos -= 1
		LN = SMILES[LNPos]
		if LN in ATOMS:
			LNindex -= 1
		if LN == ')':
			lBlockCounter += 1
		if LN == '(' and lBlockCounter != 0:
			lBlockCounter -= 1
		if LN in BONDS and lBlockCounter == 0:
			LNbond = LN
		if LN in NUMBERS:
			lOuterNumGroup = numbersHandler(LNPos)

	# Combine LN info and determine whichGroup
	if lOuterNumGroup:
		LNindicies = [LNindex, lOuterNumGroup[1], atomIndex]
		LNgroup = LN + lOuterNumGroup[0] + atomSymbol
	else:
		LNindicies = [LNindex, atomIndex]
		LNgroup = LN + atomSymbol

	info = whichGroup(LNPos, atomPosition, LNgroup, LNindicies)
	if info is not FALSE:
		FGinfo.append(info)

	# If entries are in FGinfo, return it. Otherwise return FALSE
	if FGinfo:
		return FGinfo
	else:
		return FALSE



# Creates groups stemming from CHARGES symbol within SMILEScode
def chargeHandler(chargeSymbol, chargePostion, atomIndex):

	# Charges are a special case within SMILEScode. Each individual charge is enclosed
	# in brackets with the charged atom, in the form [C-] or [N+]. Charged hygroden
	# atoms are removed from the SMILEScode by the formatSmilesCode() function, since
	# no functioanl group
	chargeGroup = SMILEScode[chargePostion-2:chargePostion+2]
	info = whichGroup(chargePostion-2, chargePostion+2, chargeGroup, [atomIndex])
	return info




# Determines the proper group from a theoretical group built by the data parsers, and returns the group info found if determined, or False if nothing was matched
def whichGroup(startPosition, endPosition, group, atomIndicies):

	# Holds templates which contain group, sent to expandGroup to attempt a fullMatch
	portionMatches = []

	# Only one fullMatch is possible
	fullMatch = FALSE

	# Loop through groups
	for line in f.open('FGlist.txt', 'r'):

		# Variables
		lineInfo = re.compile(r'\S+').findall(line)
		FGtemplate = lineInfo[0].replace('[R]', 'R')
		FGname = lineInfo[1]
		difference = len(template) - len(group)

		match = checkGroup(group, FGtemplate, FALSE)

		# If group is larger than the template, or no match was found, skip
		if difference < 0 or match is FALSE:
			continue

		# Equivalent lengths with a true match means an identical match
		elif difference == 0 and match is TRUE:
			fullMatch = TRUE
			fullDict = createFGDict(atomIndicies, FGtemplate)
			fullInfo = (FGname, FGtemplate, atomIndicies, fullDict)

		# Unequivalent lengths with a true match means a portion match
		elif difference > 0 and match is TRUE:
			portionMatches.append([FGtemplate, FGname])

	# If a group was found without any portions, return group match information
	if fullMatch and not portionMatches:
		return fullInfo

	# If a fullMatch and a portion was found, determine if on of the portions complete a group. Take expandGroup if it exists, otherwise take fullMatch
	elif fullMatch and portionMatches:
		expandedGroup = expandGroup(startPosition, endPosition, group, atomIndicies, portionMatches)
		if expandedGroup is not FALSE:
			return expandedGroup
		else:
			return fullInfo

	# If there are portions, determine if they complete a group
	elif fullMatch is FALSE and portionMatches:
		expandedGroup = expandGroup(startPosition, endPosition, group, atomIndicies, portionMatches)
		if expandedGroup is not FALSE:
			return expandedGroup
		else:
			return FALSE

	# If there is no matches at all, return FALSE
	elif fullMatch is FALSE and not portionMatches:
		return FALSE


# Expands a group passed into whichGroup if it is found to be a portion of another larger group.
# Returns fully matched and expanded group info, or False if nothing was matched
# portionMatches is a list of larger templates which contain the subtemplate, group, within it
def expandGroup(startPosition, endPosition, group, atomIndicies, portionMatches):

	for portionGroup in portionMatches:

		# Variables
		template = portionGroup[0]
		FGname = portionGroup[1]
		difference = len(template) - len(group)
		groupRegex = re.compile(group)
		posInTemplate = groupRegex.search(template).start()
		leftRequiredPos = startPosition - posInTemplate
		atomIndex = atomIndicies[0] - len(re.compile(ATOMS).findall(template[0:posInTemplate]))  # Initialize as leftmostIndex in group subtracted from the number of atoms to the left of groups appearance in the tempate
		rightRequiredPos = endPosition + len(template[(posInTemplate+len(group)):len(template)]) + 1  # Extra 1 necessary to compensate for 0 index String
		requiredGroup = SMILEScode[leftRequiredPos:rightRequiredPos]
		finalAtomIndicies = []

		numCount = len(re.compile(r'/D').findall(requiredGroup))
		if numCount != 0:
			requiredPos = -1
			formattedRequiredGroup = ""
			for symbol in requiredGroup:
				if symbol in ATOMS:
					atomIndex += 1
					if atomIndex != finalAtomIndicies[-1]: # If the index is not in finalAtomIndicies, add it
						finalAtomIndicies.append(atomIndex)
				requiredPos+=1 # This is effectively the postion of the number in requiredPos later used
				formattedRequiredGroup += symbol

				if symbol in NUMBERS:
					numInfo = numbersHandler(requiredPos)
					numGroup = '(' + numInfo[0] + ')'
					isLeft = True if requiredPos < posInTemplate else FALSE

					# numGroup must be able to "fit" inside of the fully constructed requiredGroup, where the numbers are taken into account
					if isLeft is True and (len(numGroup) - len(requiredGroup[0:requiredPos]) >= len(numGroup):
						# Get rid of number and replace it with numGroup
						formattedRequiredGroup[0:len(formattedRequiredGroup-1)] += numGroup # Add number
						finalAtomIndicies.append(numInfo[1]) # Add outer index
						formattedRequiredGroup = formattedRequiredGroup[len(numGroup):len(formattedRequiredGroup)] # Chop len(numGroup) amount of symbols to the left
						finalAtomIndicies = finalAtomIndicies[2:len(finalAtomIndicies)]

					# If isLeft is false, isRight is essentially True. Same idea, numGroup must "fit" inside constructed requiredGroup, or essentially template
					elif isLeft if FALSE and (len(numGroup) - len(requiredGroup[requiredPos:len(requiredGroup)]))
						formattedRequiredGroup[0:len(formattedRequiredGroup-1)] += numGroup
						finalAtomIndicies.append(numInfo[1]) # Add outer index
						formattedRequiredGroup += requiredGroup[requiredPos:len(requiredGroup)-len(numGroup)]
						atomIndex -=  len(re.compile(ATOMS).findall(requiredGroup[requiredPos:len(requiredGroup)]))# Account for deletion of atoms at the end, so subtract index by 2 to account for deletion od indicies

					# If a number can't be evaluted in such terms, the group expanded automatically is false
					else:
						# The number is left inside formattedRequiredGroup in this case, so match will always be FALSE if this is reached
						# because template cannot have numbers, and the positioing of the requiredGroup does not match with template
						break
			match = checkGroup(formattedRequiredGroup, template, True)
			if match is FALSE:
				continue
			else:
				finalDict = createFGDict(finalAtomIndicies, template)
				return (FGname, template, finalAtomIndicies, finalDict)

		# If no numbers are found, directly check
		elif numCount == 0:
			match = checkGroup(requiredGroup, template, True)
			if math is FALSE:
				continue
			else:
				finalAtomIndicies = []
				for symbol in requiredGroup:
					if symbol in ATOMS:
						atomIndex += 1
						finalAtomIndicies.append(atomIndex)
				finalDict = createFGDict(finalAtomIndicies, template)
				return (FGname, template, finalAtomIndicies, finalDict)
	return FALSE



# Creates external group if number is found, called only by data parsers to add to theoretical group being built
# Returns info in the form [atomSymbol, atomIndex, atomPosition]
def numbersHandler(numPosition):

	for RINGGROUP in RINGPOSITIONS:

		# If the number OPENS the group, return the CLOSING info
		if RINGGROUP[0][0] == numPosition:
			return(RINGGROUP[1][1], RINGGROUP[1][2], RINGGROUP[1][3])

		# If the number CLOSES the group, return the OPENING info
		elif RINGGROUP[1][0] == numPosition:
			return(RINGGROUP[0][1], RINGGROUP[0][2], RINGGROUP[0][3])


# Creates a FG dictionary which correlates maingroup and rgroup atoms to their proper indicies within the strucutre
def createFGDict(FGindicies, FGtemplate):
	FGdict = []
	indexCounter = -1
	for symbol in FGtemplate:
		if symbol in ATOMS or symbol == 'R':
			indexCounter += 1
			FGdict.append([symbol, FGindicies[indexCounter]])
	return FGdict


# Returns TRUE if a template and a group are the same, or if a group is contained within the template
# Always returns FALSE if len(group) > len(template)
# if FULL=TRUE, only matches group === template, but not group within template
# if FULL=FALSE, matches group within template as well
def checkGroup(group, template, full=TRUE):

	difference = len(template) - len(group)

	# Group larger than template, return FALSE
	if difference < 0:
		return FALSE

	# If only fullMatch is possible, but difference exists, return FALSE
	if difference > 0 and full is TRUE:
		return FALSE

	# Group same length as template
	if difference == 0:
		groupCounter = -1
		for symbol in template:
			groupCounter += 1
			groupSymbol = group[groupCounter]
			if groupSymbol != symbol:
				return FALSE
		return TRUE

	# Group smaller than template
	if difference > 0 and full is FALSE:
		for shift in range(0, difference):
			templateCounter = - 1 + shift
			symbolCounter = - 1
			for symbol in group:
				symbolCounter += 1
				templateCounter += 1
				templateSymbol = template[templateCounter]
				if templateSymbol != symbol:
					break
				if symbolCounter == len(group) - 1:
					return TRUE
		return FALSE


# Only works for single charged hydrogens
def formatSmilesCode(smiles):
	smilesPos = -1
	for symbol in smiles:
		smilesPos += 1
		if symbol == '[':
			startBracketPos = smilesPos
		if symbol == 'H':
			cutPos = smilesPos
			while smiles[cutPos] != '+':
				cutPos += 1
			reFormatted = smiles[0:startBracketPos-1] + smiles[startBracketPos+1] + smiles[cutPos+1:len(smiles)]
			break
