""" Use this file to test the SMILES decoder on a test SMILES code to create the software molecule object
    See its decoded data placed into data structures
"""

from pandas import read_excel
from molecule import Molecule

def commonKeys(a,b):
    common_keys_set = set()
    for i in a.keys():
            for j in b.keys(): 
                if i == j:
                    if not i in common_keys_set:
                        common_keys_set.add(i)
    return common_keys_set

def genDiffList(): 

    ifg_original_set = read_excel("output/ref.xlsx", sheet_name="Precise Functional Groups")
    ifg_original_set = ifg_original_set.set_index("Refcode")

    diff_list: "list[tuple]" = []
    for line in open('resources/smiles.txt'):
        _, smiles, name = [x for x in line.strip().split(" ") if x]
        mol = Molecule(smiles, name, "mol")
        d = ifg_original_set.loc[name]
        ref_fgs = {fg: fg_count for (fg, fg_count) in dict(d[2:-5]).items() if fg_count != 0}
        mol_fgs = mol.functional_groups_exact

        if ref_fgs != mol_fgs:
            diff_list.append((smiles, name, mol_fgs, ref_fgs))
    return diff_list


def testSingleMol(smiles, name):
    
    ifg_original_set = read_excel("output/ref.xlsx", sheet_name="All Functional Groups")
    ifg_original_set = ifg_original_set.set_index("Refcode")

    d = ifg_original_set.loc[name]
    ref_fgs = {fg: fg_count for (fg, fg_count) in dict(d[2:-5]).items() if fg_count != 0}
    # print(ifg_original_set.columns[2:-5])
    mol = Molecule(smiles, name, "mol")
    mol_fgs = mol.functional_groups_all
    print("mol_fgs = ", mol_fgs)
    print("ref_fgs = ", ref_fgs)
    

diff_list = genDiffList()
for v in diff_list:
    print("smiles = ", v[0])
    print("name = ", v[1])
    print("mol_fgs = ", v[2])
    print("ref_fgs = ", v[3])
print(len(diff_list))


# testSingleMol('COC(=O)C1CNc2ccccc2N2C=CN=C2C1O', 'APENIT')

# # # mol_test = Molecule("OCC1OC(CC1O)N1C=NC2=C1C(=O)NC=N2", "APENIT", "mol")
# mol_test = Molecule("CC(=C)C1CCC23OC2C(CC2(C)CC(=O)C(CC(=O)C1)O2)OC3=O", "YEGVUC", "mol")
# print(mol_test.ring_atoms)
# print(mol_test.edges)
# # print(mol_test.size)
# # for edge in mol_test.edges:
# #     print(edge.index)
# print(mol_test.functional_groups_all)
# print(mol_test.functional_groups_exact)
# fg_test = Molecule("RC(OR)OR", "HemiAcetal", "fg")
# for v in fg_test.verticies:
#     print(v)
#     print("explicit = ", v.explicit_degree)
#     print("implicit = ",v.implicit_degree)
#     print("total = ", v.total_degree)
# # from datetime import datetime
# # now_t = datetime.now()
# diff_list: "list[str]" = []


# symbol = 'Br+'
# charge = '+'
# print(''.join([char for char in symbol if charge != char]))

# for i in range(0,10):
#     if i == 5:
#         i = 20
#     print(i)

# from molecule import Molecule
# from datetime import datetime
# now_t = datetime.now()
# for line in open('resources/smiles.txt'):
#     _, smiles, name = [x for x in line.strip().split(" ") if x]
#     mol = Molecule(smiles, name, "mol")
# print(datetime.now() - now_t)
#     # d = ifg_original_set.loc[name]
#     # ref_fgs = {fg: fg_count for (fg, fg_count) in dict(d[2:-5]).items() if fg_count != 0}
#     # mol_fgs = mol.functional_groups_all
    # if ref_fgs != mol_fgs:
    #     diff_list.append(name)

# ifg_original_set = read_excel("output/ref.xlsx", sheet_name="Precise Functional Groups")
# ifg_revised_set = read_excel("output/FunctionalGroups.xlsx", sheet_name="Precise Functional Groups")
# print(ifg_original_set.equals(ifg_revised_set))
# diff = concat([ifg_original_set,ifg_revised_set]).drop_duplicates(keep=False)
# if not diff.empty: print(diff)
    # print(d)

    # #### Ring Validation #####
    # if name in ["BUDWIH", "EPOMAY", "GILXUV", "MAMHEO", "ODEJOX", "OHIXUY", "OYEXIA", "ROCLUR", "TOLNAK"]:
    #     continue
    # print("Aromatic: ", mol.aromatic_ring_count,  "==", d["aromaticRingCount"])
    # print("Non-Aromatic: ", mol.non_aromatic_ring_count,  "==", d["nonAromaticRingCount"])
    # print("Total:", mol.total_rings, "==", d["aromaticRingCount"] + d["nonAromaticRingCount"])
    # assert mol.aromatic_ring_count == d["aromaticRingCount"]
    # assert mol.non_aromatic_ring_count == d["nonAromaticRingCount"]
    # assert mol.total_rings == d["aromaticRingCount"] + d["nonAromaticRingCount"]
    # print('\nBy passing object of class')
    # print(mol.__dict__)
# print(diff_list)
# mol = Molecule("RC(R)(O)OR", "AFAFAP", "fg")
# # # "O=C1CCC2=C(C#N)c3cc4OCOc4cc3CCN12"
# print(mol.__dict__)
# print("verticies:")
# print(mol.verticies)
# print("")
# print("order:")
# print(mol.order)
# print("")
# print("edges:")
# print([str(e) for e  in mol.edges])
# print("")
# print("size: ")
# print(mol.size)
# print("")
# print("atom_counts: ")
# print(mol.atom_counts)

