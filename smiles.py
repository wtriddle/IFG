import ifg
import re

# for line in open('SMILES.txt', 'r'):
#     lineInfo = re.compile(r'\S+').findall(line)
#     smiles = lineInfo[0]
#     print("EVALUATING ", lineInfo[1])
#     functionalGroupData = ifg.ifg(smiles)
#     print(functionalGroupData)
functionalGroupData = ifg.ifg("COC(=O)C([n+]1ccccc1)=C(C(=O)OC)[C-]1C(=O)COC1=O")
print(functionalGroupData)
print("COC(=O)C([n+]1ccccc1)=C(C(=O)OC)[C-]1C(=O)COC1=O")

list2 = [15,16,17,18,19]
list1 = [16,17,18]
print("list1 = ", list1)
print("list2 = ", list2)
result1 = all(index in list1 for index in list2)
result2 = all(index in list2 for index in list1)
print("result1 = ", result1)
print("result2 = ", result2)
