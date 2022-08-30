"""Constant Variables used across the program source files"""

import re

##### Regular Expressions #####
ATOM_REGEX = re.compile(r'Br[+-]?|Cl[+-]?|[a-zA-Z][+-]?')
CHARGE_REGEX = re.compile(r'[+-]{1}')
BOND_REGEX = re.compile(r'[=#]')
DIGIT_REGEX = re.compile(r'[0-9]{1}')
PARENTH_REGEX = re.compile(r'[()]')
BRACKET_REGEX = re.compile(r'[\[\]]')
SMILES_REGEX = re.compile(r'(Br[+-]?|Cl[+-]?|[a-zA-Z][+-]?|[=#]|[0-9]{1}|[()])')
AMINO_ACID_REGEX = re.compile(r'[nN]H[23]?\+')

##### Valence Electrons Per-Atomic Orbital Counts #####
VALENCE_COUNTS = {
    'C': 4,
    'N': 3,
    'P': 3,
    'O': 2,
    'S': 2,
    'SE': 2,
    'F': 1,
    'CL': 1,
    'BR': 1,
    'I': 1,
    'R': 1
}

##### Electrons Per-Bond Counts #####
ELECTRON_BOND_COUNTS = {
    '': 1,
    '=': 2,
    '#': 3
}