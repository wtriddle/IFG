# IFG (Identify Functional Groups)

Python algorithms which extract functional group data from SMILES (Simplified molecular-input line-entry system) codes <br>
Project and Research Overview: https://youtu.be/feT36AVNRgk
<p align="center">
  <!-- <img src="https://github.com/wtriddle/IFG/blob/master/MoleculeGifSmall.gif" />
  <img src="https://github.com/wtriddle/IFG/blob/master/MoleculeGifSmall.gif" /> -->
</p>

# Installation

A detailed install video for beginners is here:
https://www.youtube.com/watch?v=kU3W452HC9s

Otherwise, follow these steps: <br>

1. Download anaconda and git <br>
https://www.anaconda.com/products/individual <br>
https://git-scm.com/downloads <br>

2. Download the repository to your computer with:
```
  git clone https://github.com/wtriddle/IFG.git
```
3. Set your PYTHONPATH to include the /src folder <br>
On Windows, edit environemnt variables and create PYTHONPATH to your /src location <br>
On Linux, open .bashrc and place the following line
```
  export PYTHONPATH="Your/Path/To/Ifg/Source:$PYTHONPATH"
```
4. Create a folder called output in the root directory
5. Activate the base anaconda environment in the root directory of IFG with
```
  conda activate base
```
6. Run index.py in the anaconda base environemnt to test the program

# Usage

Index.py is a top-level script which handles the src files, but those files can be directly used as well. This algorithm provides a decoded digital model of molecular SMILES codes. The implementation of IFG itself is an extension of this digital molecular model.

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
    print(functionalGroups.allFgs)
    print(functionalGroups.preciseFgs)
```

# Configuration
A detailed configuration video for beginners is here: pending <br>

Otherwise, follow the instructions below

## Adding Functional Groups
Go into the FGlist.txt file and add the template with the name of the functional group to the end of the text file. <br>
Ensure that the format contains [R] for R group atoms and that the atoms included in the new FG are all presente in Molecule <br>

## Choosing a new set of SMILES codes
Go into <b> config.py <b> and edit the ``` SMILESPATH ``` variable to target the new SMILES text list. <br>
Ensure that the new set of SMILES codes is disjointed and follows the form for each line
```
  SMILES REFOCDE
```

# Contributing

When contributing, open an issue for suggestions and follow the commenting style observed across the repository <br>
New formats of SMILES codes can be supported by updating and contributing to ``` Molecule ``` <br>
Optimization of the DFS algorithm can be changed in ``` ifg ```

# License
[MIT](https://choosealicense.com/licenses/mit/)
