import ifg

for smiles in open('SMILES.txt', 'r'):
    functionalGroupData = ifg.ifg(smiles)
    print(functionalGroupData)
