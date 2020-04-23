from molecule import molecule
from ifg import ifg
import re
import time

# elapsedtime = time.process_time()
# for line in open('smiles.txt', 'r'):
# 	lineInfo = re.compile(r'\S+').findall(line)
# 	smiles = lineInfo[1]
# 	refcode = lineInfo[2]
# 	mol = molecule(smiles)
# 	print(mol)
# 	print(refcode)
# 	print(smiles)
# elapsedtime = time.process_time()
# print(elapsedtime)
smiles = ifg("COC(=O)C=CC1C2CCC3C2C(C)(C)CCCC13C","BIVLEY")
# print(smiles.functionalGroups)
for group in smiles.functionalGroups:
	print(group.NAME, group.SMILES, group.atomData)
# firstFG = smiles.functionalGroups[0]
# print(firstFG.SMILES)
# print(firstFG.atomData)
# template = molecule("RC(=O)OR")

# for i in range(0,len(template.atomData)):
# 	print(i, template.atomData[i], template.bondData[i])

# print('\n')

# for tempAtom in template.atomData:
# 	index = tempAtom[0]
# 	if tempAtom[1] == 'R':
# 		template.bondData[index].clear()

# for i in range(0,len(template.atomData)):
# 	print(i, template.atomData[i], template.bondData[i])

# smiles.removeBond(1,template)
# template.atomData[1][0] = 2

# print('\n')

# for i in range(0,len(template.atomData)):
# 	print(i, template.atomData[i], template.bondData[i])
# print('\n')
# for i in range(0,len(smiles.atomData)):
# 	print(i, smiles.atomData[i], smiles.bondData[i])
# expandGroup = smiles.expandGroup(smiles.atomData[2],1,template)
# print(expandGroup)
# print(template.atomData)
# print("I Atom 	      Bond")
# print('\n')
# for i in range(0,len(template.atomData)):
# 	print(i, template.atomData[i], template.bondData[i])

# for atom in template.atomData:
# 	index = atom[0]
# 	if atom[1] == 'R':
# 		template.bondData[index].clear()

# print('\n')
# for i in range(0,len(template.atomData)):
# 	print(i, template.atomData[i], template.bondData[i])

# atomIndex = 1
# for index,bonds in enumerate(template.bondData):
# 	for j,bond in enumerate(bonds):
# 		if bond[0] == atomIndex:
# 			del(template.bondData[index][j])
# print('\n')
# for i in range(0, len(template.atomData)):
# 	print(i, molSmi.atomData[i], molSmi.bondData[i])

# print('\n')
# for i in range(0,len(template.atomData)):
# 	print(i, template.atomData[i], template.bondData[i])
