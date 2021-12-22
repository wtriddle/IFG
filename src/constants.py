import re

# Symbol matching lists and regex's for decoding process
ATOM_REGEX = re.compile(r'[a-zA-Z]')
CHARGE_REGEX = re.compile(r'\+|\-')
BOND_REGEX = re.compile(r'\=|\#')

# Double Lettered atom to Single Letter Atom representation
DLA_TO_SAR = {                                  
    "Br": "X",                                       
    "Cl": "Z"
}

# Charecter Matching Lists
ATOMS =['C', 'O', 'N', 'S', 'I','F', 'P',           # SAR atoms non-aromatic
        'c', 'n', 'o','s', 'i', 'f', 'p',           # SAR atoms aromatic
        'R',                                        # R group atom (used for FG templates)
        'X', 'Z',                                   # DLA->SAR converted atoms
        ]
BONDS = ['=', '#']
BRACKETS = ['[', ']']
CHARGES = ['+', '-']
NUMBERS = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
NON_BRANCHING_SYMBOLS = ATOMS + BONDS + BRACKETS + CHARGES