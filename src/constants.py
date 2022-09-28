"""Constant variables used across the program source files"""

import re

##### Regular Expressions #####
ATOM_REGEX = re.compile(r'Br[+-]?|Cl[+-]?|[a-zA-Z][+-]?')
"""Regular Expression that captures all atom symbols in a hydrogen-suppressed SMILES code, inclusive of charge and double lettered atoms"""
CHARGE_REGEX = re.compile(r'[+-]{1}')
"""Regular Expression that captures a charged symbol"""
BOND_REGEX = re.compile(r'[=#]')
"""Regular Expression that captures a bond symbol"""
DIGIT_REGEX = re.compile(r'[0-9]{1}')
"""Regular Expression that captures a digit symbol"""
PARENTH_REGEX = re.compile(r'[()]')
"""Regular Expression that captures a parenthetical symbol"""
BRACKET_REGEX = re.compile(r'[\[\]]')
"""Regular Expression that captures a bracketed symbol"""
SMILES_REGEX = re.compile(r'(Br[+-]?|Cl[+-]?|[a-zA-Z][+-]?|[=#]|[0-9]{1}|[()])')
"""Regular Expression that captures all individual symbols of significance in a hydrogen-suppressed SMILES code"""
AMINO_ACID_REGEX = re.compile(r'[nN]H[23]?\+')
"""Regular Expression that matches an amino acid"""

##### Required Valence Electrons Per-Atomic Orbital To Fulfill Electron Configuration Counts #####
REQUIRED_VALENCE_COUNTS: "dict[str, int]" = {
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
"""Dictionary of atomic symbol to valence electrons required"""

##### Electrons Per-Bond Counts #####
ELECTRON_BOND_COUNTS: "dict[str, int]" = {
    '': 1,
    '=': 2,
    '#': 3
}
"""Dictionary of bond symbol to number of valence electrons given"""