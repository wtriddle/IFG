# Written by William Riddle from 9/23/19 -

import re
import collections
import sys

# The data parsing functions written here slice the SMILEScode in various ways to obtain a string of a potential functionalGroup, represented in the FGlist.txt
# To optimze the program, the he data parsing functions can find the same group multiple times, such that a single functionalgroup may appear multiple times in FGdata
# This is intended so that the SMILEScode can be looked at from various viewpoints, using its linearity

# VARIBLE DESCRIPTIONS

# RN stands for rightNeighbor, and evalutes symbols to the right of the symbol of interest in the SMILEScode code
# LN stands for leftNeighbor, and evalutes symbols to the left of the symbol of interest in the SMILEScode code

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

# GITposition = global variable containing the start position of a group within a template
# i.e. C(=O) in RC(=O)R would be 1
# Data appended in checkGroup, and resets for every individual symbol.
# Relevant in expandGroup

RINGPOSITIONS = []
FGdata = []
GITposition = []

# Constants
ATOMS = ['C','O','N','c','n','o']
RGROUP = ['C','O','N','c','n','o']
CHARGES = ['+','-']
PARENTHESIS = ['(', ')']
BRACKETS = ['[', ']']
BONDS = ['=', '#']
NUMBERS = ['1','2','3','4','5','6','7','8','9']
ATOMSREGEX = re.compile(r'[a-zA-Z]')
NUMBERSREGEX = re.compile(r'[0-9]+')
RGROUPREGEX = re.compile(r'R')

# Head function to loop through SMILEScode and call proper data parsers on given characters. Final FGdata is formatted as well
def ifg(SMILES):
	global SMILEScode
	SMILEScode = formatSMILEScode(SMILES) # Definition of SMILEScode
	# #print("EVALUATING ", SMILEScode)
	global SMILEScodelength
	SMILEScodelength = len(SMILEScode)
	SMILEScodePos = -1
	atomIndex = 0
	initializeRINGPOSITIONS()
	for symbol in SMILEScode:
		# reset GITposition each new symbol
		GITposition.clear()
		SMILEScodePos += 1

		if symbol in PARENTHESIS:
			# #print("SMILEScode[", SMILEScodePos, "] = ", symbol, " and belongs to the PARENTHESIS group")
			continue

		elif symbol in CHARGES:
			#print("SMILEScode[", SMILEScodePos, "] = ", symbol, " and belongs to the CHARGE group")
			chargeGroup = chargeHandler(symbol, SMILEScodePos, atomIndex)
			#print("chargeHandler found ", chargeGroup, "\n")
			if chargeGroup is not False:
				FGdata.append(chargeGroup)
			del(chargeGroup)

		elif symbol in ATOMS and SMILEScodePos != 0 and SMILEScodePos != SMILEScodelength - 1:
			atomIndex += 1
			#print("SMILEScode[", SMILEScodePos, "] = ", symbol, " and belongs to the ATOMS group, and is atom number ", atomIndex)
			directGroup = atomHandler(symbol, SMILEScodePos, atomIndex)
			#print("atomHandler found ", directGroup, "\n")
			if directGroup is not False:
				for group in directGroup:
					FGdata.append(group)
			del(directGroup)

		elif symbol in BONDS:
			#print("SMILEScode[", SMILEScodePos, "] = ", symbol, " and belongs to the BONDS group")
			bondGroup = bondHandler(symbol, SMILEScodePos, atomIndex)
			#print("bondHandler found ", bondGroup, "\n")
			if bondGroup is not False:
				FGdata.append(bondGroup)
			del(bondGroup)

	FGdataFinal = evaluateFGdata()
	#print(FGdataFinal)
	FGdataDict = {}
	for group in FGdataFinal:
		keys = []
		# #print("Group is ", group)
		# #print(group[0])
		for key in FGdataDict.keys():
			keys.append(key)
		# #print(keys)
		if group[0] not in keys:
			# #print(group[0], " is not in ", keys)
			FGdataDict.update({group[0] : 1})
		else:
			FGdataDict[group[0]] += 1
		del(keys)
	#print("Functional Groups Final Result:\n\n")
	#print(FGdataDict)
	GITposition.clear()
	FGdata.clear()
	RINGPOSITIONS.clear()
	return (FGdataFinal, FGdataDict)


# Initializes the global RINGPOSITIONS according to its model
def initializeRINGPOSITIONS():

	# Initialize variables
	SMILEScodePos = atomIndex = -1
	evaluatedNumbers = []
	chargeGroup = ""

	# Main loop
	for symbol in SMILEScode:
		#print(RINGPOSITIONS)
		SMILEScodePos += 1

		# Keep track of ATOM positional data for when a number is found
		if symbol in ATOMS:
			atomIndex += 1
			# Atom found before number is opening atom of the ring.
			openingAtom = symbol
			openingAtomPos = SMILEScodePos

		# When an unevalutaed number is found, evalute it
		if symbol in NUMBERS and SMILEScodePos not in evaluatedNumbers:
			#print("Found ", symbol, " at position ", SMILEScodePos)
			#print(evaluatedNumbers)
			# SMILEScodePos pointing at the position of the atom symbol in SMILEScode 1 space before the first occurance of the number

			# However, if it is not an atom in this scope, then it must be a bracket surrounding a charged atom.
			# The following if statement below is a special case where charged atom is opening the ring.
			# A closing bracket must precede the position of the number.
			if SMILEScode[SMILEScodePos-1] == ']' and SMILEScodePos >= openingAtomPos + 2:
				specialCounter = SMILEScodePos - 1 # Initialized at posion of bracket
				specialSymbol = SMILEScode[specialCounter]
				while specialSymbol != '[':
					chargeGroup += specialSymbol
					specialCounter -= 1
					specialSymbol = SMILEScode[specialCounter]
				chargeGroup += '['
				chargeGroup = chargeGroup[::-1]
				# Reassign openingAtom and openingAtomPos so that the data can be \
				# appended to RINGPOSITIONS with a single statement
				openingAtom = chargeGroup
				openingAtomPos = specialCounter

			# subLoop counters
			correlatingPos = SMILEScodePos
			correlatingIndex = atomIndex

			# subLoop to find correlating atom in ring
			for subSymbol in SMILEScode[SMILEScodePos+1:len(SMILEScode)]:
				correlatingPos += 1

				# Keep track of ATOM positional data for when the correlating atom is found
				if subSymbol in ATOMS:
					correlatingAtom = subSymbol
					correlatingIndex += 1
					correlatingAtomPos = correlatingPos

				# If the same number which began the sub loop is found, then append both opening and closing info
				if subSymbol == symbol:
					#print("Correlating number ", subSymbol, " was found at position ", correlatingPos)

					# If the ring opens with a charge, it must close with a charge. Re-assign positional variables
					if chargeGroup != "":
						correlatingChargeGroup = ""
						specialCounter = correlatingPos -1
						specialSymbol = SMILEScode[specialCounter]
						while specialSymbol != '[':
							correlatingChargeGroup += specialSymbol
							specialCounter -= 1
							specialSymbol = SMILEScode[specialCounter]
						correlatingChargeGroup = correlatingChargeGroup[::-1]
						correlatingAtom = correlatingChargeGroup
						correlatingAtomPos = specialCounter

					RINGPOSITIONS.append([
						[
						SMILEScodePos,
						openingAtom,
						atomIndex,
						openingAtomPos],
						[correlatingPos,
						correlatingAtom,
						correlatingIndex,
						correlatingAtomPos]
						])

					# Add numbers to evaluted list to prevent reevalutaion of correlatingAtom in main loop
					evaluatedNumbers.append(SMILEScodePos)
					evaluatedNumbers.append(correlatingPos)
					break

				# Error caught if a correlatingPos runs out of the scope of the SMILEScode
				if correlatingPos == len(SMILEScode) - 1:
					print("There has been an error in the creation of the RINGPOSITIONS global used in the program\n")
					print("The SMILEScode code ", SMILEScode, " had an error")
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
# If no comparitive conclusion is made, then group is an individual group found by the data parsers.
# D. group will be added if:
# (1) compareGroup above was not added
# (2) group is not already within formattedFGdata
def evaluateFGdata():
	formattedFGdata = []
	groupCounter = -1
	#print("FG data is ", FGdata)
	addGroup = True
	for group in FGdata:
		addGroup = True
		#print("Evalutating", group)

		groupCounter += 1

		# group variables
		groupTemplate = group[1]
		groupIndicies = group[2]
		groupDict = group[3]
		compareCounter = -1

		# Check if group is already in formattedFGdata (A.1) Skip if True
		inFormattedData = isGroupInList(formattedFGdata, group)
		if inFormattedData is True:
			continue

		for compareGroup in FGdata:

			# compareGroup variables
			compareCounter += 1
			compareTemplate = compareGroup[1]
			compareIndicies = compareGroup[2]
			compareDict = compareGroup[3]

			# Check if loop is comparing the exact same group (A.2) Skip if True
			if compareCounter == groupCounter:
				continue

			# Determine number of Maingroup atoms
			numMainGroupAtoms = len(groupIndicies) - len(RGROUPREGEX.findall(groupTemplate))
			numMainCompareAtoms = len(compareIndicies) - len(RGROUPREGEX.findall(compareTemplate))

			# Determine Maingroup crossover
			mainGroupCrossover = 0
			# #print(groupDict)
			# #print(compareDict)
			for atom in groupDict:
				for compareAtom in compareDict:
					# #print(atom, atom[0], atom[1])
					# #print(compareAtom, compareAtom[0], compareAtom[1])
					if atom[0] == compareAtom[0] and atom[0] != 'R' and compareAtom[0] != 'R': # If the symbols are symbol equivalent and both main group Atoms
						if int(atom[1]) == int(compareAtom[1]): # If the atoms are the same index
							mainGroupCrossover += 1 # Then atoms must be of the exact same maingroup in the structure

			#print("mainGroupCrossover between ", groupTemplate, " and ", compareTemplate, " is ", mainGroupCrossover)
			#print("numMainGroupAtoms = ", numMainGroupAtoms)
			#print("numMainCompareAtoms = ", numMainCompareAtoms)
			# Check if compareGroup is fully contained within group (B.1) If ture, skip this compareGroup.
			if mainGroupCrossover == numMainCompareAtoms and len(groupIndicies) > len(compareIndicies):
				continue

			# Check if group is fully contained within compareGroup (C.1) If so, and is a new group, add compareCounter
			if mainGroupCrossover == numMainGroupAtoms and len(groupIndicies) < len(compareIndicies):
				addGroup = False
				inFormattedData = isGroupInList(formattedFGdata, compareGroup)
				if inFormattedData is False:
					formattedFGdata.append(compareGroup)
				break
		# If the group is not already in formattedFGdata, and it is not contained within another group (D.1, D.2) Add the group
		if addGroup is True:
			formattedFGdata.append(group)
		#print(formattedFGdata)

	return formattedFGdata


# Called from evaluateFGdata, checks if a group is already contained within the list.
# Only EXACT matches, i.e fginfo === fginfo for all cells is True
def isGroupInList(formattedFGdata, group):

	# If there is an atom in formattedFGdata
	if formattedFGdata:
		# #print("formattedFGdata is ", formattedFGdata)

		for compareGroup in formattedFGdata:
			# #print("Compare group is ", compareGroup)

			# Comparative Variables
			compareGroupTemplate = compareGroup[1]
			compareGroupIndicies = compareGroup[2]
			compareGroupDict = compareGroup[3]
			numMainCompareAtoms = len(compareGroupIndicies) - len(RGROUPREGEX.findall(compareGroupTemplate))

			# Group Variables
			groupTemplate = group[1]
			groupIndicies = group[2]
			groupDict = group[3]
			numMainGroupAtoms = len(groupIndicies) - len(RGROUPREGEX.findall(groupTemplate))

			# Check for equivalent number of atoms
			if len(compareGroupIndicies) == len(groupIndicies):
				mainGroupCrossover = 0
				for atom in groupDict:
					for compareAtom in compareGroupDict:

						# Check if atoms are maingroup equivalent, and if they are in the same positional index.
						# The atoms in both functionalGroup scopes are "identical" if so, so add to mainGroupCrossover
						if atom[0] == compareAtom[0] and atom[0] != 'R':
							if int(atom[1]) == int(compareAtom[1]):
								mainGroupCrossover += 1

				# If there is full maingroups crossover and the templates are the same
				if mainGroupCrossover == numMainCompareAtoms == numMainGroupAtoms and compareGroupTemplate == groupTemplate:
						# Then the exact group is already contained within formattedFGdata
						return True
			else: # if not equivalent lengths, continue
				continue
	return False


# Creates groups stemming from BONDS symbol within SMILEScode
def bondHandler(bondSymbol, bondPostion, atomIndex):

	# Initialize positional variables to create group
	LNPos = bondPostion - 1
	RNPos = bondPostion + 1
	LNindex = RNindex = atomIndex
	LN = SMILEScode[LNPos]
	RN = SMILEScode[RNPos]


	# Special case which occurs commonly as R(=R)
	if LN == '(' and SMILEScode[LNPos-1] in ATOMS and RN in ATOMS and SMILEScode[RNPos+1] == ')':
		atomIndicies = [atomIndex, atomIndex+1]
		bondGroup = SMILEScode[LNPos-1] + "(=" + RN + ')'
		# #print("Bond group from special case is ", bondGroup)
		FGinfo = whichGroup(LNPos-1, RNPos+1, bondGroup, atomIndicies)
		if FGinfo is not False:
			return FGinfo
		else:
			return False


	# LN loop
	lBlockCounter = 0
	lOuterParenthGroup = ""
	lOuterNumGroup = ""


	# Find the first non blocking group atom connected to the bond
	# If LN is an atom, loop not run
	while LN not in ATOMS and lBlockCounter > 0:
		if LN == ')':
			lBlockCounter += 1
		if LN == '(' and lBlockCounter != 0:
			lBlockCounter -= 1
		if LN in NUMBERS and lBlockCounter == 0:
			lOuterNumGroup = numbersHandler(LNPos)
		if lBlockCounter != 0:
			lOuterParenthGroup += symbol
		if LN in ATOMS:
			LNindex -= 1
		if lBlockCounter != 0 and LN in ATOMS:
			lOuterIndicies.append(LNindex)
		LNPos -= 1
		LN = SMILEScode[LNPos]
	else:
		LNindex -= 1

	# RN loop. If RN is already an atom, loop not run
	# The atom next to the bond symbol is always the atom involved in the bond, no blocking groups occur
	while RN not in ATOMS:
		if RN in ATOMS:
			RNindex += 1
		RNPos += 1
		RN = SMILEScode[RNPos]
	else:
		RNindex += 1

	# lOuterParenthGroup evaluation to see if the group should be included in bondGroup
	lOuterParenthGroup = lOuterParenthGroup[::-1] # Flip because appended from the left, so group appears backwards
	atomCount = len(ATOMSREGEX.findall(lOuterParenthGroup))
	numCount = len(NUMBERSREGEX.findall(lOuterParenthGroup))
	if atomCount > 1 or numCount != 0:
		lOuterParenthGroup = ""

	# Bond Group creation according to the specific case
	if lOuterNumGroup and lOuterParenthGroup:
		atomIndicies = [LNindex, lOuterNumGroup[1], lOuterIndicies[0], RNindex]
		bondGroup = LN + lOuterNumGroup[0] + lOuterParenthGroup + bondSymbol + RN
	elif lOuterNumGroup and not lOuterParenthGroup:
		atomIndicies = [LNindex, lOuterNumGroup[1], RNindex]
		bondGroup = LN + lOuterNumGroup[0] + bondSymbol + RN
	elif lOuterParenthGroup and not lOuterNumGroup:
		atomIndicies = [LNindex, lOuterIndicies[0], RNindex]
		bondGroup = LN + lOuterParenthGroup + bondSymbol + RN
	else:
		atomIndicies = [LNindex, RNindex]
		bondGroup = LN + bondSymbol + RN

	# Determine the group and return the info
	FGinfo = whichGroup(LNPos, RNPos, bondGroup, atomIndicies)
	return FGinfo



# Creates groups stemming from ATOMS symbol within SMILEScode
# Special RN case for first inner atom + outer atom group not written yet. Place in if needed.
def atomHandler(atomSymbol, atomPosition, atomIndex):
	# #print("atomHandler is creating initial functional group shells...\n")
	# Initialize Variables
	LNPos = atomPosition - 1
	RNPos = atomPosition + 1
	RN = SMILEScode[RNPos]
	LN = SMILEScode[LNPos]
	RNindex = LNindex = atomIndex
	RNbond = LNbond = ""
	FGinfo = []

	# ----RN LOOP----
	# #print("-----RN LOOP-----")
	# #print("INTIAL CONDITIONS:\n")
	# #print("rBlockCounter = ", rBlockCounter)
	# #print("rInnerParenthGroup = ", rInnerParenthGroup, "\n")
	# #print("SMILEScode[", RNPos, "] = RN = ", RN)
	# #print("RNindex = ", atomIndex)
	rBlockCounter = 0
	rOuterGroups = 0
	rInnerParenthGroup = "("
	rInnerParenthIndicies = []
	rInnerParenthGroupPositions = []
	numGroup = ""

	# Loop continues until an atom within the same grouping layer as atomSymbol is found to the right
	loopCounter = 0
	rOuterBoolean = True
	while rBlockCounter > 0 or RN not in ATOMS:
		#print(rOuterGroups)
		if RNPos == len(SMILEScode) - 1:
			RN = ""
			#print("NO Right group exists...")
			break
		loopCounter += 1
		if RN in ATOMS:
			RNindex += 1
		# Add RN numbers handler
		# #print(loopCounter)
		# IF RN in the same grouping layer as atomSymbol

		# If scope of RN loop is within a parenthesis, capture the information inside of it with respect to atomSymbol
		if rOuterBoolean is False and rBlockCounter == 0:
			#print("Reset rOuterBoolean to ", rOuterBoolean)
			rOuterBoolean = True
		if rBlockCounter != 0 and rOuterBoolean is True:
			rOuterGroups += 1
			rOuterBoolean = False
			#print("Found inner group and added one to ", rOuterGroups, " and set rOuterBoolean to ", rOuterBoolean)
		if rBlockCounter != 0:
			rInnerParenthGroup += RN
			rInnerParenthGroupPositions.append(RNPos)
		if rBlockCounter != 0 and RN in ATOMS:
			rInnerParenthIndicies.append(RNindex)
		if rBlockCounter != 0 and RN in NUMBERS:
			numGroup = numbersHandler(RNPos)
			rInnerParenthIndicies.append(numGroup[1])
			rInnerParenthGroup += numGroup[0]
		#
		# #print("RN LOOP INFO:\n")
		# #print("rBlockCounter = ", rBlockCounter)
		# #print("RNpos = ", RNPos)
		# #print("SMILEScode[", RNPos, "] = RN = ", RN)
		# #print("rInnerParenthGroup = ", rInnerParenthGroup)
		# #print("RNindex = ", RNindex)
		# #print("numGroup = ", numGroup)
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
		# #print("Made it to the end of this loop!")
		# #print("RN = ", RN)
		# #print("InnerLoop = ", innerLoop)
	else:
		RNindex += 1
		#print("CEASED ON ", loopCounter)
		#print("\nRN LOOP RESULT:\n")
		#print("rBlockCounter = ", rBlockCounter)
		#print("RNpos = ", RNPos)
		#print("SMILEScode[", RNPos, "] = RN = ", RN)
		#print("rInnerParenthGroup = ", rInnerParenthGroup)
		#print("RNindex = ", RNindex)
		#print("numGroup = ", numGroup)
		#print("Not terminated by a break!")
	rInnerAtoms = ATOMSREGEX.findall(rInnerParenthGroup)

	# #print("Matches = ", matches)
	atomCount = len(rInnerAtoms)

	#print("atomCount == ", atomCount)
	#print(rOuterGroups)
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
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[-1], rInnerParenthGroup, [atomIndex, rInnerParenthIndicies[0]])
		if info is not False:
			FGinfo.append(info)
		del(info)
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[0], rInnerParenthGroup, [atomIndex, rInnerParenthGroup[0], RNindex])
		if info is not False:
			FGinfo.append(info)
	elif atomCount == 1 and rOuterGroups < 1:

		rInnerParenthGroup = atomSymbol + rInnerParenthGroup
		# #print("Evalutating rInnerParenthGroup", rInnerParenthGroup)
		info = whichGroup(atomPosition, rInnerParenthGroupPositions[-1]+1, rInnerParenthGroup, [atomIndex, rInnerParenthIndicies[0]])
		if info is not False:
			FGinfo.append(info)
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
		firstParenthesisGroupIndicies = rInnerParenthIndicies[0:firstEndAtom]
		secondParenthesisGroup = '(' + rInnerParenthGroup[firstEndParenthesis:len(rInnerParenthGroup)]
		secondParenthesisGroupIndicies = rInnerParenthIndicies[firstEndAtom:len(rInnerParenthIndicies)]
		#print(firstParenthesisGroup, " with indicies ", firstParenthesisGroupIndicies)
		#print(secondParenthesisGroup, " with indicies ", secondParenthesisGroupIndicies)
		if firstParenthesisGroup != secondParenthesisGroup:
			rInnerParenthGroupV1 = atomSymbol + '(' + firstParenthesisGroup[1] + ')' + secondParenthesisGroup[1]
			info = whichGroup(atomPosition, rInnerSmilesPos + 1, rInnerParenthGroupV1, [atomIndex, firstParenthesisGroupIndicies[0], secondParenthesisGroupIndicies[0]])
			if info is not False:
				FGinfo.append(info)
			rInnerParenthGroupV2 = atomSymbol + '(' + secondParenthesisGroup[1] + ')' + secondParenthesisGroup[1]
			info = whichGroup(atomPosition, rInnerParenthGroupPositions[0], rInnerParenthGroupV2, [atomIndex, secondParenthesisGroup[0], firstParenthesisGroup[0]])
			if info is not False:
				FGinfo.append(info)
	elif rOuterGroups > 2:
		#print("There was an error")
		return 0

	# Combine RN info and determine whichGroup. If RN is empty, do not run this portion of code
	# #print("RN = ", RN)
	if RN != "":
		RNgroup = atomSymbol + RNbond + RN
		print("RNgroup is ", RNgroup)
		RNindicies = [atomIndex, RNindex]
		info = whichGroup(atomPosition, RNPos, RNgroup, RNindicies)
		if info is not False:
			FGinfo.append(info)
		print("found, ", info)


	# ----LN LOOP----
	lOuterNumGroup = ""
	lBlockCounter = 0

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
		if LN in NUMBERS:
			lOuterNumGroup = numbersHandler(LNPos)
		LNPos -= 1
		LN = SMILEScode[LNPos]
	else:
		LNindex -= 1

	# Combine LN info and determine whichGroup
	if lOuterNumGroup:
		LNindicies = [LNindex, lOuterNumGroup[1], atomIndex]
		LNgroup = LN + LNbond + '(' +  lOuterNumGroup[0] + ')' + atomSymbol
		print("LNgroup with number is ", LNgroup)
	else:
		LNindicies = [LNindex, atomIndex]
		LNgroup = LN + LNbond + atomSymbol


	info = whichGroup(LNPos, atomPosition, LNgroup, LNindicies)
	if info is not False:
		FGinfo.append(info)

	# #print("atomHandler found the group info to be, ", FGinfo)
	# If entries are in FGinfo, return it. Otherwise return False
	if FGinfo:
		return FGinfo
	else:
		return False



# Creates groups stemming from CHARGES symbol within SMILEScode
def chargeHandler(chargeSymbol, chargePostion, atomIndex):

	# Charges are a special case within SMILEScode. Each individual charge is enclosed
	# in brackets with the charged atom, in the form [C-] or [N+]. Charged hygroden
	# atoms are removed from the SMILEScode by the formatSMILEScode() function, since
	# no functioanl group
	chargeGroup = SMILEScode[chargePostion-2:chargePostion+2]
	info = whichGroup(chargePostion-2, chargePostion+2, chargeGroup, [atomIndex])
	return info




# Determines the proper group from a theoretical group built by the data parsers, and returns the group info found if determined, or False if nothing was matched
def whichGroup(startPosition, endPosition, group, atomIndicies):


	# #print("finding the group for ", group)
	# Holds templates which contain group, sent to expandGroup to attempt a fullMatch
	portionMatches = []

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
		# 	#print("match between ", group, " and ", FGtemplate, " is, ", match)
		# If group is larger than the template, or no match was found, skip
		if difference < 0 or match is False:
			continue

		# Equivalent lengths with a True match means an identical match
		elif difference == 0 and match is True:
			fullMatch = True
			fullDict = createFGDict(atomIndicies, FGtemplate)
			fullInfo = (FGname, FGtemplate, atomIndicies, fullDict)

		# Unequivalent lengths with a True match means a portion match
		elif difference > 0 and match is True:
			portionMatches.append([FGtemplate, FGname])

	# If a group was found without any portions, return group match information
	if fullMatch and not portionMatches:
		print("fullMatch True, portionMatches False")
		return fullInfo

	# If a fullMatch and a portion was found, determine if on of the portions complete a group. Take expandGroup if it exists, otherwise take fullMatch
	elif fullMatch and portionMatches:
		expandedGroup = expandGroup(startPosition, endPosition, group, atomIndicies, portionMatches)
		if expandedGroup is not False:
			print("Fullmatch with Expanded Group")
			return expandedGroup
		else:
			print("Fullmatch with NO Expand Group")
			return fullInfo

	# If there are portions, determine if they complete a group
	elif fullMatch is False and portionMatches:
		expandedGroup = expandGroup(startPosition, endPosition, group, atomIndicies, portionMatches)
		if expandedGroup is not False:
			print("No full match but Expand Group")
			return expandedGroup
		else:
			print("No full match or expandGroup")
			return False

	# If there is no matches at all, return False
	elif fullMatch is False and not portionMatches:
		print("No full match or portionMatches")
		return False


# Expands a group passed into whichGroup, if it is found to be a portion of another larger group.
# Returns fully matched and expanded group info, or False if nothing was matched
# portionMatches is a list of larger templates which contain the subtemplate, group, within it
# startPosition = string index of first symbol of a group within SMILEScode
# endPosition =
def expandGroup(startPosition, endPosition, group, atomIndicies, portionMatches):

	# #print("Expand group called with ", portionMatches)
	portionCounter = -1
	for portionGroup in portionMatches:
		portionCounter += 1
		# Variables
		template = portionGroup[0]
		FGname = portionGroup[1]
		difference = len(template) - len(group)
		posInTemplate = GITposition[portionCounter]
		# #print(len(GITposition), " != ", len(portionMatches))
		if len(GITposition) != len(portionMatches):
				# #print("GITposition error! ")
				return 0
		# #print("GITposition is ", GITposition)
		# #print("position of ", group," in ", template, " is ", posInTemplate)
		leftRequiredPos = startPosition - posInTemplate
		# #print("leftRequiredPos is ", leftRequiredPos)
		atomIndex = atomIndicies[0]
		# #print("template is being cut from ", intialPosition, " to ", finalPosition)
		groupToRight = template[posInTemplate + len(group):len(template)] # Symbols to the right of group within template.
		# #print("Group to right is ", groupToRight)
		rightRequiredPos = endPosition + len(groupToRight) + 1 # Add one for zero index
		leftRequiredGroup = SMILEScode[leftRequiredPos:startPosition]
		numLeftAtoms = len(ATOMSREGEX.findall(leftRequiredGroup))
		rightRequiredGroup = SMILEScode[endPosition+1:rightRequiredPos]
		numRightAtoms = len(ATOMSREGEX.findall(rightRequiredGroup))
		# #print(leftRequiredPos, " is compared to ", 0)
		# #print(rightRequiredPos, " is compared to ", len(SMILEScode) - 1)
		if (leftRequiredPos < 0) or (rightRequiredPos > (len(SMILEScode) - 1)): # If the group is outside the scope of the SMILEScode, then group cannot possible exist
			# #print("SKIPPED")
			continue
		# requiredGroup works for noncontinguous SMILEScode groups
		requiredGroup = leftRequiredGroup + group + rightRequiredGroup
		finalAtomIndicies = []
		#print("\nVARIABLES FOR THE PORTION GROUP: \n", portionGroup)
		#print("atomIndicies = ", atomIndicies)
		#print("template = ", template)
		#print("group = ", group)
		#print("leftRequiredGroup = ", leftRequiredGroup)
		#print("rightRequiredGroup = ", rightRequiredGroup)
		#print("requiredGroup = ", requiredGroup)
		#print("difference = ", difference)
		#print("posInTemplate = ", posInTemplate)
		#print("startPosition = ", startPosition)
		#print("rightRequiredPos = ", rightRequiredPos)
		#print("leftRequiredPos = ", leftRequiredPos)
		#print("endPosition = ", endPosition)
		#print("numLeftAtoms = ", numLeftAtoms)
		#print("numRightAtoms = ", numRightAtoms)
		#print("Right required template = ", template[posInTemplate+len(group):len(template)])
		#print("Expanding ", group, " into ", portionGroup, "where in the smiles code, the template is ", requiredGroup)
		numCount = len(re.compile(r'/D').findall(requiredGroup))
		# #print("numCount = ", numCount)
		if numCount != 0:
			requiredPos = -1
			formattedRequiredGroup = ""
			for symbol in requiredGroup:
				if len(formattedRequiredGroup) == len(template):
					break
				if symbol in ATOMS:
					atomIndex += 1
					if atomIndex != finalAtomIndicies[-1]: # If the index is not in finalAtomIndicies, add it
						finalAtomIndicies.append(atomIndex)
				requiredPos+=1 # This is effectively the postion of the number in requiredGroup later used
				formattedRequiredGroup += symbol
				# #print("formattedRequiredGroup is ", formattedRequiredGroup)

				if symbol in NUMBERS:
					numInfo = numbersHandler(requiredPos)
					numGroup = '(' + numInfo[0] + ')'
					isLeft = True if requiredPos < posInTemplate else False

					# numGroup must be able to "fit" inside of the fully constructed requiredGroup, where the numbers are taken into account
					if isLeft is True and ((len(requiredGroup[0:requiredPos+1]) - len(numGroup)) >= 0):
						# Get rid of number and replace it with numGroup
						formattedRequiredGroup = formattedRequiredGroup[0:len(formattedRequiredGroup)-1] + numGroup # Add number
						finalAtomIndicies.append(numInfo[1]) # Add outer index
						numAtomsToChop = len(ATOMSREGEX.findall(formattedRequiredGroup[0:len(numGroup)]))
						formattedRequiredGroup = formattedRequiredGroup[len(numGroup):len(formattedRequiredGroup)] # Chop len(numGroup) amount of symbols from the left
						finalAtomIndicies = finalAtomIndicies[numAtomsToChop:len(finalAtomIndicies)]

					# If isLeft is False, isRight is essentially True. Same idea, numGroup must "fit" inside constructed requiredGroup, or essentially template
					elif isLeft is False and (len(requiredGroup[requiredPos:len(requiredGroup)]) - len(numGroup) >= 0):
						formattedRequiredGroup = formattedRequiredGroup[0:len(formattedRequiredGroup)-1] + numGroup
						finalAtomIndicies.append(numInfo[1]) # Add outer index
					# If a number can't be evaluted in such terms, the group expanded automatically is False
					else:
						# The number is left inside formattedRequiredGroup in this case, so match will always be False if this is reached
						# because template cannot have numbers, and the positioing of the requiredGroup does not match with template
						break
			# #print("EXPANDED ", group, " INTO ", formattedRequiredGroup)
			if len(formattedRequiredGroup) == len(template):
				match = checkGroup(formattedRequiredGroup, template, True)
			else:
				print("Expand Group ERROR")
				sys.exit()
			# #print("MATCH IS ", match)
			if match is False:
				continue
			else:
				finalDict = createFGDict(finalAtomIndicies, template)
				GITposition.clear()
				return (FGname, template, finalAtomIndicies, finalDict)

		# If no numbers are found, directly check
		elif numCount == 0:
			# #print("Checking ", requiredGroup)
			match = checkGroup(requiredGroup, template, True)
			# #print(requiredGroup, " is ", match, " with ", template)
			if match is False:
				continue
			else:
				finalAtomIndicies = []
				#print("atomIndicies = ", atomIndicies)
				#print("numLeftAtoms = ", numLeftAtoms)
				#print("numRightAtoms = ", numRightAtoms)
				for number in range(atomIndicies[0], atomIndicies[0] - numLeftAtoms - 1, -1):
					finalAtomIndicies.insert(0, number)
				#print(finalAtomIndicies)
				for index in atomIndicies[1:len(atomIndicies)-1]:
					finalAtomIndicies.append(index)
				#print(finalAtomIndicies)
				for number in range(atomIndicies[-1], atomIndicies[-1] + numRightAtoms + 1):
					# #print(finalAtomIndicies)
					finalAtomIndicies.append(number)
				#print(finalAtomIndicies)
					# #print(finalAtomIndicies)
				# #print("Final atom indicies is ", finalAtomIndicies)
				# #print(atomIndicies[-1] + numRightAtoms)
				# #print(atomIndicies[-1])
				finalDict = createFGDict(finalAtomIndicies, template)
				GITposition.clear()
				return (FGname, template, finalAtomIndicies, finalDict)
	GITposition.clear()
	return False



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


# Returns True if a template and a group are the same, or if a group is contained within the template
# Always returns False if len(group) > len(template)
# if FULL=True, only matches group === template, but not group within template
# if FULL=False, matches group within template as well
def checkGroup(group, template, full=True):

	group = group.upper()
	difference = len(template) - len(group)
	# #print("Checking ", group, " aginast ", template)
	# #print("There is a difference of ", difference, " while full = " , full)
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


# Only works for single charged hydrogens
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
			break
	if reFormatted == "":
		return SMILEScode
	else:
		return reFormatted
