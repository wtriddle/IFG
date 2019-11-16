import ifg
import re

for line in open('SMILES.txt', 'r'):
    lineInfo = re.compile(r'\S+').findall(line)
    smiles = lineInfo[0]
    functionalGroupData = ifg.ifg(smiles)
    print(functionalGroupData)
    break
