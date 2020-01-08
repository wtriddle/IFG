import ifg
import re
from openpyxl import Workbook

# Openpyxl objects
# fgWorkbook = Workbook()
# fgSheet = fgWorkbook.active

# Set intial column names
# fgSheet.cell(row=1, column=1).value = "Refcode"
# fgSheet.cell(row=1, column=2).value = "SMILES"

# columnCounter = 2
# rowCounter = 1
# fullDictionary = {'SecondaryAmine': 376, 'Ketone': 424, 'Alkyne': 30, 'TetiaryAmine': 210,
# 'Nitroso': 200, 'Ether': 423, 'Alcohol': 365, 'Ester': 274, 'PrimaryAmine': 69,
# 'Carboxylate': 16, 'Aldehyde': 55, 'Nitrile': 64, 'SecondaryKetimine': 150,
# 'Alkene': 461, 'Amide': 70, 'Carbamate': 10, 'SecondaryAldimine': 74, 'Imide': 14,
# 'Azide': 19, 'PrimaryKetimine': 14, 'Isocynate': 3, 'Azo': 23, 'PrimaryAldimine': 21,
# 'Hydroperoxide': 5, 'Isonitrile': 2, 'Nitrosooxy': 8, 'Carbonate': 3, 'Peroxide': 2}

# for group in fullDictionary.items():
#     groupName = group[0]
#     columnCounter += 1
#     fgSheet.cell(row=1, column=columnCounter).value = groupName

# for line in open('smiles.txt', 'r'):
#     rowCounter += 1
    # for column in range(1, columnCounter+1):
    #     fgSheet.cell(row=rowCounter, column=column).value = 0
    # lineInfo = re.compile(r'\S+').findall(line)
    # smilesWH = lineInfo[0]
    # smilesWOH = lineInfo[1]
    # RefCode = lineInfo[2]
    # print(lineInfo)
    # fgSheet.cell(row=rowCounter, column=1).value = RefCode
    # fgSheet.cell(row=rowCounter, column=2).value = smiles
    # print("EVALUATING ", lineInfo[0], " ", smilesWOH)
    # functionalGroupData = ifg.ifg(smilesWOH, smilesWH)
    # print(functionalGroupData)
    # for group in functionalGroupData[1].items():
    #     for column in range(1,columnCounter+1):
    #         if fgSheet.cell(row=1, column=column).value == group[0]:
    #             fgSheet.cell(row=rowCounter, column=column).value = int(group[1])
    #             break
# fgWorkbook.save("functionalGroupData.xlsx")
# print(fullDictionary)
#
data = ifg.ifg("COc1ccc(cc1)C1=NC(C)(OC21C(=O)N(C)c1ccccc21)C(=O)OC", "[CH3]Oc1[cH][cH]c([cH][cH]1)C1=NC([CH3])(OC21C(=O)N([CH3])c1[cH][cH][cH][cH]c21)C(=O)O[CH3]")
print(data)
