import ifg
import re

# for line in open('SMILES.txt', 'r'):
#     lineInfo = re.compile(r'\S+').findall(line)
#     smiles = lineInfo[0]
#     print("EVALUATING ", lineInfo[1])
#     functionalGroupData = ifg.ifg(smiles)
#     print(functionalGroupData)
functionalGroupData = ifg.ifg("CC12CCC(O)C3(C)C1C(OC2=O)C=C1COC(=O)C=C31")
print(functionalGroupData)
print("CC12CCC(O)C3(C)C1C(OC2=O)C=C1COC(=O)C=C31")
