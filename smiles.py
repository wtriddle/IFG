import ifg
import re
from openpyxl import Workbook

# Openpyxl objects
# fgWorkbook = Workbook()
# fgSheet = fgWorkbook.active
#
# # Set intial column names
# fgSheet.cell(row=1, column=1).value = "Refcode"
# fgSheet.cell(row=1, column=2).value = "SMILES"
#
#
# # List which holds functional group names
# functionalGroups = []
#
# # Grab all names from FGlist.txt
# for line in open('FGlist.txt', 'r'):
#     lineInfo = re.compile(r'\S+').findall(line)
#     groupName = lineInfo[1]
#     if groupName not in functionalGroups:
#         functionalGroups.append(groupName)
#
# functionalGroups.append("Alcohol")
# functionalGroups.append("Acetal")
# functionalGroups.append("Hemiketal")
# functionalGroups.append("Hemiacetal")
#
# # Sort them alphabetically
# functionalGroups.sort()
# holder = functionalGroups.copy()
#
# # Add cyclic/aromatic distinctions to list
# indexCounter = insertions = 0
# for name in holder:
#     indexCounter += 1
#     functionalGroups.insert(indexCounter+insertions, "Cyclic" + name)
#     if len(re.compile(r'Amine').findall(name)) != 0 or name == "Alcohol":
#         functionalGroups.insert(indexCounter+insertions, "Aromatic" + name)
#         insertions += 1
#     insertions += 1
#
# # Insert names into the top row
# for column in range(3, len(functionalGroups)+3):
#     fgSheet.cell(row=1,column=column).value = functionalGroups[column-3]
#
# fgSheet.cell(row=1, column=len(functionalGroups)+3).value = "aromaticRingCount"
# fgSheet.cell(row=1, column=len(functionalGroups)+4).value = "nonAromaticRingCount"
# fgSheet.cell(row=1, column=len(functionalGroups)+5).value = "RingCount"
# fgSheet.cell(row=1, column=len(functionalGroups)+6).value = "AminoAcid"
#
# rowCounter = 1
# maxColumn = fgSheet.max_column
# print(functionalGroups)
#
# for line in open('smiles.txt', 'r'):
#     rowCounter += 1
#     for column in range(2, maxColumn+1):
#         fgSheet.cell(row=rowCounter, column=column).value = 0
#     lineInfo = re.compile(r'\S+').findall(line)
#     smiles = lineInfo[1]
#     RefCode = lineInfo[2]
#     fgSheet.cell(row=rowCounter, column=1).value = RefCode
#     fgSheet.cell(row=rowCounter, column=2).value = smiles
#     print("EVALUATING ", lineInfo[2], " ", smiles)
#     functionalGroupData = ifg.ifg(smiles)
#     print(functionalGroupData)
#     for group in functionalGroupData[1].items():
#          for column in range(2,maxColumn+1):
#             if fgSheet.cell(row=1, column=column).value == group[0]:
#                 fgSheet.cell(row=rowCounter, column=column).value = int(group[1])
#                 if group[0] == "AminoAcid" and int(group[1]) == 1:
#                     fgSheet.cell(row=rowCounter, column=column).value = "Yes"
#                 elif group[0] == "AminoAcid" and int(group[1]) == 0:
#                     fgSheet.cell(row=rowCounter, column=column).value = "No"
#                 break
#     del(functionalGroupData)
# fgWorkbook.save("functionalGroupData.xlsx")

data = ifg.ifg("CC1(C)OC2OC(C(O)C2O1)C(=O)C1CC2CC1C=C2")
for thing in data:
    print("\n")
    print(thing)
print("CC1(C)OC2OC(C(O)C2O1)C(=O)C1CC2CC1C=C2")

# zero = 0
# print(zero)
# zero = ['entry']
# print(zero)
# myString = "SecondaryAmine"
# myRegex = re.compile(r'(?<=Secondary)\S+')
# match = myRegex.search(myString)
# print(myString)
# print(myRegex)
# print(match)
# print(match.group(0))

# file = open("FGheirarchy.txt", "r+")
#
#
# for line in open('FGheirarchy.txt', 'r'):
#     lineInfo = re.compile(r'\S+').findall(line)
#     heirarchy = lineInfo[0].split(":")
#     print(heirarchy)
#     if "PrimaryAmine" in heirarchy:
#         print("Yessitmate")
#     else:
#         print("snowsir!nosnowsinthesnowsir!")
