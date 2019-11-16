import ifg
import re

for line in open('SMILES.txt', 'r'):
    lineInfo = re.compile(r'\S+').findall(line)
    smiles = lineInfo[0]
    # posCounter = -1
    # print(" The length of the smiles code is ", len(smiles))
    # for symbol in smiles:
    #     posCounter += 1
    #     print(symbol, " is in position ", posCounter)
    #     print("Chopped smiles string is ", smiles[0:posCounter+1])
    functionalGroupData = ifg.ifg(smiles)
    print(functionalGroupData)
    break
