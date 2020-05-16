# IFG (Identify Functional Groups)
Python scripts which extract functional group data from SMILES (Simplified molecular-input line-entry system) codes

molecule.py creates object molecules based on the SMILES code representation, which contain bonding information about each atom, ring information, and other specifications such as the precense of an amino acid or alcohol

ifg.py uses the molecule object from molecule.py to determine functional groups based on a CHON smiles code 

