# IFG (Identify Functional Groups)

A python package which identifies functional groups in
canonicalized non-isomeric non-isotopic hydrogen-suppressed
[SMILES](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html) codes 
(The Simplified Molecular-Input Line-Entry System) <br>

This package has been documented at the [IFG docs](https://wtriddle.github.io/IFG/) online through Github Pages. Please view this documentation 
to get started using this package for identification of functional groups in SMILES codes on your computer.

To get a quick start on the IFG python package and its research goals, please view the
[Project and Research Overview Video](https://youtu.be/yOdvyQ0seAc).

# Installation

To install IFG for SMILES-based functional group data collection and data analysis, follow the steps below:

1. Install [git](https://git-scm.com/downloads)
2. Install [poetry](https://python-poetry.org/docs/)
3. Clone this github repository into a desired folder with:

   ```cmd
     C:\> git clone https://github.com/wtriddle/IFG.git
   ```

4. Install the python packages required for IFG with poetry:

   ```cmd
     C:\IFG> poetry install
   ```

5. Add the absolute IFG/ifg path to your system PYTHONPATH environemnt variable (example C:\IFG\ifg).
   See this [tutorial](https://www.simplilearn.com/tutorials/python-tutorial/python-path) to set this up on your computer

You are now ready to use the IFG package! Please read the
[IFG docs](https://wtriddle.github.io/IFG/) to view how to use the python package in your own work.

# Examples

The python code examples below show the SMILES-based functional group identification in action using the IFG python package defined in this repository:

```python
    >>> from chem.molecule import Molecule
    >>> mol = Molecule('O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O', name='ABEGOH', type="mol")
    >>> print(mol)
    O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O
    >>> mol.functional_groups_all
    {
      'Non Aromatic Ketone': 1,
      'Non Aromatic Amide': 2,
      'Non Aromatic SecondaryAmine': 1,
      'Non Aromatic TertiaryAmine': 3,
      'Nitro': 3
    }
    >>> mol.non_aromatic_ring_count
    2
```

```python
    >>> from chem.molecule import Molecule
    >>> mol = Molecule('CC(=O)Oc1ccccc1', name="AXUDIH", type="mol")
    >>> print(mol)
    CC(=O)Oc1ccccc1
    >>> mol.functional_groups_all
    {'Ketone': 1, 'Ester': 1, 'Ether': 1}
    >>> mol.aromatic_ring_count
    1
```

# Contributing

If you would like to contribute to the repository by including support of different SMILES codes, or by incorporating new SMILES-based computation, or by including new features, please open an issue in this repository about your topic. It will be seen by the maintainer of this project and reviewed, then allowed or rejected.

If an error or problem with the software arises, please open an issue so the maintainer of this project can review it and understand the problem.
Be as detailed as possible and make sure it is relevant to the code which has already been developed. If you have developed a plugin or
included a new feature that you would like to see work with the code, please open an issue and still inform the maintainer of this.

# License

[MIT](https://choosealicense.com/licenses/mit/)


# Doc Relaunch
Planning on relaunching docs through github pages again.
