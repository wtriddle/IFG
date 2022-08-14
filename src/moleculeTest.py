""" Use this file to test the SMILES decoder on a test SMILES code to create the software molecule object
    See its decoded data placed into data structures
"""


from Molecule import Molecule
for line in open('resources/smiles.txt'):
    _, smiles, name = [x for x in line.strip().split(" ") if x]
    print(smiles, name)
    mol = Molecule(smiles, name)
    print('\nBy passing object of class')
    print(mol.__dict__)                                 # __dict__ is the winner
"""
The only atoms apart of a NPR are those which appear linealry after it has been opnened as the MROR and the atoms which open and close the ring
the maximum number of polycyclic atoms for a NPR is 2, and these appear at the end of 

MROR (NPR/PR) (OPEN/CLOSE) (NPR/PR)
0 00
0 01
0 10
0 11
1 00
1 01
1 10
1 11 
OPEN/CLOSE
NPR PR
8 different cases/behaviors


"""
# mol = Molecule("COc1ccc(cc1)C1CC(=O)CC(=C(C)C)C(C1)=C(C)C", "AFAFAP")
# "O=C1CCC2=C(C#N)c3cc4OCOc4cc3CCN12"
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
# print("rings: ")
# print(mol.rings)
# print("")
# print("ring_counts: ")
# print(mol.ring_counts)
# print("")
# print("atom_counts: ")
# print(mol.atom_counts)

