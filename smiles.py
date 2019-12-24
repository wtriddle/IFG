import ifg
import re

for line in open('SMILES.txt', 'r'):
    lineInfo = re.compile(r'\S+').findall(line)
    smiles = lineInfo[0]
    print("EVALUATING ", lineInfo[1], " ", smiles)
    functionalGroupData = ifg.ifg(smiles)
    print(functionalGroupData)
# functionalGroupData = ifg.ifg("N#Cc1cccc(Oc2cccnc2)c1C#N")
# print(functionalGroupData)
# print("N#Cc1cccc(Oc2cccnc2)c1C#N")
