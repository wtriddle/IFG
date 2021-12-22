""" Use this file to test the SMILES decoder on a test SMILES code to create the software molecule object
    See its decoded data placed into data structures
"""


from Molecule import Molecule

mol = Molecule("[NH3+]CCC(O)C(=O)[O-]", "RABTUL")
print('\nBy passing object of class')
print(mol.__dict__)                                 # __dict__ is the winner

print("atomData:")
print(mol.atomData)
print("")
print("bondData:")
print(mol.bondData)
print("")
print("AROMATICINDICES: ")
print(mol.AROMATICINDICES)
print("")
print("CYCLICINDICES: ")
print(mol.CYCLICINDICES)
print("")
print("RING_CLOSE_POSITIONS: ")
print(mol.RING_CLOSE_POSITIONS)
print("")
print("RING_COMPLEMENTS: ")
print(mol.RING_COMPLEMENTS)
print("")
print("RING_SELF: ")
print(mol.RING_SELF)
print("")
print("RING_OPEN_POSITIONS: ")
print(mol.RING_OPEN_POSITIONS)
print("")
print("ringData: ")
print(mol.ringData)
print("")
print("ALCOHOLICINDICES: ")
print(mol.ALCOHOLICINDICES)
print("")

