from ifg import ifg
from ifgTest import ifgTest
import re
import sys
from helpers import createFgDataDict
import os

for line in open(os.getcwd() + '/src/resources/smiles.txt', 'r'):
    lineInfo = re.compile(r'\S+').findall(line)
    smiles = lineInfo[1]
    refcode = lineInfo[2]
    print(smiles)
    if smiles == "C#CN(c1ccccc1)c1ccc(cc1)N(C#C)c1ccccc1" or smiles == "CCCCN1c2ccccc2c2ccc3cc(OC)ccc3c12":
        continue
    base = ifg(smiles, refcode)
    test = ifgTest(smiles, refcode)
    baseAllDict = createFgDataDict(base.functionalGroups)
    basePreciseDict = createFgDataDict(base.preciseFunctionalGroups)
    testAllDict = createFgDataDict(test.functionalGroups)
    testPreciseDict = createFgDataDict(test.preciseFunctionalGroups)
    print("Base: ", baseAllDict)
    print("Altered: ", testAllDict)
    if baseAllDict != testAllDict:
        print("all fgs not equal")
        sys.exit()

    if basePreciseDict != testPreciseDict:
        print("precise fgs not equal")
        sys.exit()
