# IFG (Identify Functional Groups)

Python algorithms which extract functional group data from SMILES (Simplified molecular-input line-entry system) codes
<p align="center">
  <!-- <img src="https://github.com/wtriddle/IFG/blob/master/MoleculeGifSmall.gif" />
  <img src="https://github.com/wtriddle/IFG/blob/master/MoleculeGifSmall.gif" /> -->
</p>

# Installation

To generate the organic functional groups from a molecule, the following packages are required:

```
    xlsxwriter
    pandas
    progress
```

To execute the program, run python version 3.5+ in an anaconda environemnt with the above packages on the index.py file.
An excel file of the generated smiles codes, according those given in the smiles.txt under resources, will be output and pre-cleaned, thanks
to pandas

# Usage

This algorithm provides a decoded digital model of molecular SMILES codes. The implementation of IFG itself is an extension of this digital molecular model.

Digital Molecule:
```
    from Molecule import Molecule
    mol = Molecule('O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O', 'ABEGOH')
    print(mol)
    >> ABEGOH : O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O
```


Identification of functional groups in a SMILES code:
```
    from ifg import ifg
    functionalGroups = ifg(SMILES="NN1C=NN=C1N",REFCODE="VUPTAC02")
```

# Contributing

When contributing, open an issue for suggestions and follow the commenting style observed across the repository

# License
[MIT](https://choosealicense.com/licenses/mit/)
