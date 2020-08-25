import ifg
import ifgTest as testing
import re
import sys


def createFgDataDict(functionalGroups):
    fgDataDict = {}
    for group in functionalGroups:
        if group.NAME in fgDataDict.keys():
            fgDataDict[group.NAME] += 1
        else:
            fgDataDict.update({group.NAME: 1})
    return fgDataDict


for line in open('./../src/resources/smiles.txt', 'r'):
    lineInfo = re.compile(r'\S+').findall(line)
    smiles = lineInfo[1]
    refcode = lineInfo[2]
    base = ifg.ifg(smiles, refcode)
    test = testing.ifgTest(smiles, refcode)
    baseAllDict = createFgDataDict(base.functionalGroups)
    basePreciseDict = createFgDataDict(base.preciseFunctionalGroups)
    testAllDict = createFgDataDict(test.functionalGroups)
    testPreciseDict = createFgDataDict(test.preciseFunctionalGroups)
    if baseAllDict != testAllDict:
        print("all fgs not equal")
        print(test.functionalGroups)
        print(base.functionalGroups)
        sys.exit()

    if basePreciseDict != testPreciseDict:
        print("precise fgs not equal")
        sys.exit()
