from Molecule import Molecule
from MoleculeTest import MoleculeTest
from ifgTest import ifgTest

mol = MoleculeTest("CCCCN1c2ccccc2c2ccc3cc(OC)ccc3c12", "Test")
fgs = ifgTest("CCCCN1c2ccccc2c2ccc3cc(OC)ccc3c12", "test")
for group in fgs.functionalGroups:
    print(group.getSymbolDict())
print("CCCCN1c2ccccc2c2ccc3cc(OC)ccc3c12")
print(mol.AROMATICINDICES)
print(mol.CYCLICINDICES)

# import re
# import time
# import sys

# start = time.time()
# for line in open('../src/resources/smiles.txt', 'r'):

#     lineInfo = re.compile(r'\S+').findall(line)
#     smiles = lineInfo[1]
#     refcode = lineInfo[2]

#     mol = Molecule(smiles, refcode)
#     molTest = MoleculeTest(smiles, refcode)
#     if not mol.AROMATICINDICES == molTest.AROMATICINDICES:
#         print("Test failed 1")
#         print(mol.SMILES)
#         sys.exit()
#     if not mol.CYCLICINDICES == molTest.CYCLICINDICES:
#         print("Test failed 2")
#         print(mol.SMILES)
#         sys.exit()
#     if not mol.ringCount == molTest.ringCount:
#         print()
#         print("Test failed 3")
#         print(mol.SMILES)
#         sys.exit()
#     if not mol.aromaticCount == molTest.aromaticCount:
#         print("Test failed 4")
#         print(mol.SMILES)
#         sys.exit()
#     if not mol.nonAromaticCount == molTest.nonAromaticCount:
#         print("Test failed 5")
#         print(mol.SMILES)
#         sys.exit()
#     if not mol.RINGDICT == molTest.RINGDICT:
#         print("Test failed 6")
#         print(mol.SMILES)
#         sys.exit()
#     if not mol.atomData == molTest.atomData:
#         print("Test failed 7")
#         print(molTest.atomData)
#         print(mol.atomData)
#         print(mol.SMILES)
#         sys.exit()
#     if not mol.bondData == molTest.bondData:
#         print("Test failed 8")
#         print(mol.SMILES)
#         sys.exit()


# print(time.time() - start)
# print("Test Succeeded")

# for z, x in enumerate([1, 2, 3, 4, 5, 6, 6, 7]):
#     print(x)
#     x = 15
#     print(x)
#     if z == 0:
#         x = 15
