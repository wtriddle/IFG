import ifg
import re

for line in open('SMILES.txt', 'r'):
    lineInfo = re.compile(r'\S+').findall(line)
    smiles = lineInfo[0]
    print("EVALUATING ", lineInfo[1])
    functionalGroupData = ifg.ifg(smiles)
    print(functionalGroupData)
# functionalGroupData = ifg.ifg("CN1C(=O)N(C)c2nc3cccc(C)c3nc2C1=O")
# print(functionalGroupData)
