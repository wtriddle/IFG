# import pandas as pd
# pd.set_option('display.max_rows', 1000)

# df1 = pd.read_excel('./../src/FgTesting.xlsx')
# df2 = pd.read_excel('./../output/FunctionalGroups.xlsx')
# difference = df1[df1 != df2]
# print(difference.count())

from Molecule import Molecule

mol = Molecule(
    "[N-]=[N+]=NC(N=[N+]=[N-])c1cc(C(N=[N+]=[N-])N=[N+]=[N-])c(cc1C(N=[N+]=[N-])N=[N+]=[N-])C(N=[N+]=[N-])N=[N+]=[N-]", "Test")
for atom in mol.atomData.items():
    print(atom)
print(mol.atomData)
