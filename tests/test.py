import pandas as pd
pd.set_option('display.max_rows', 1000)

df1 = pd.read_excel('./../src/FgTesting.xlsx').sort_index(inplace=True)
df2 = pd.read_excel(
    './../output/FunctionalGroups.xlsx').sort_index(inplace=True)
difference = df1[df1 != df2]
print(difference.count())

# from Molecule import Molecule

# mol = Molecule(
#     "[N-]=[N+]=NC(N=[N+]=[N-])c1cc(C(N=[N+]=[N-])N=[N+]=[N-])c(cc1C(N=[N+]=[N-])N=[N+]=[N-])C(N=[N+]=[N-])N=[N+]=[N-]", "Test")
# for atom in mol.atomData.values():
#     print(atom)
# for atom in mol.bondData.values():
#     print(atom)
# print(mol.atomData)
