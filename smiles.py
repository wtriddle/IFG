import ifg
import re
from openpyxl import Workbook

# Openpyxl objects
fgWorkbook = Workbook()
fgSheet = fgWorkbook.active

# Set intial column names
fgSheet.cell(row=1, column=1).value = "Refcode"
fgSheet.cell(row=1, column=2).value = "SMILES"


# List which holds functional group names
functionalGroups = []

# Grab all names from FGlist.txt
for line in open('FGlist.txt', 'r'):
    lineInfo = re.compile(r'\S+').findall(line)
    groupName = lineInfo[1]
    if groupName not in functionalGroups:
        functionalGroups.append(groupName)

# Sort them alphabetically
functionalGroups.sort()
holder = functionalGroups.copy()

# Add cyclic/aromatic distinctions to list
indexCounter = insertions = 0
for name in holder:
    indexCounter += 1
    functionalGroups.insert(indexCounter+insertions, "Cyclic" + name)
    if len(re.compile(r'Amine').findall(name)) != 0:
        functionalGroups.insert(indexCounter+insertions, "Aromatic" + name)
        insertions += 1
    insertions += 1

# Insert names into the top row
for column in range(3, len(functionalGroups)):
    fgSheet.cell(row=1,column=column).value = functionalGroups[column-3]


rowCounter = 1
maxColumn = fgSheet.max_column
# print(functionalGroups)

for line in open('smiles.txt', 'r'):
    rowCounter += 1
    for column in range(2, maxColumn+1):
        fgSheet.cell(row=rowCounter, column=column).value = 0
    lineInfo = re.compile(r'\S+').findall(line)
    smilesWH = lineInfo[0]
    smilesWOH = lineInfo[1]
    RefCode = lineInfo[2]
    fgSheet.cell(row=rowCounter, column=1).value = RefCode
    fgSheet.cell(row=rowCounter, column=2).value = smilesWOH
    # print("EVALUATING ", lineInfo[0], " ", smilesWOH)
    functionalGroupData = ifg.ifg(smilesWOH, smilesWH)
    # print(functionalGroupData)
    for group in functionalGroupData[1].items():
        for column in range(2,maxColumn+1):
            if fgSheet.cell(row=1, column=column).value == group[0]:
                fgSheet.cell(row=rowCounter, column=column).value = int(group[1])
                break
fgWorkbook.save("functionalGroupData.xlsx")
