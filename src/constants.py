import re

ATOM_REGEX = re.compile(r'Br[+-]?|Cl[+-]?|[a-zA-Z][+-]?')
CHARGE_REGEX = re.compile(r'\+|\-')
BOND_REGEX = re.compile(r'[=#]')
DIGIT_REGEX = re.compile(r'[0-9]{1}')
PARENTH_REGEX = re.compile(r'[()]')
BRACKET_REGEX = re.compile(r'[\[\]]')
SMILES_REGEX = re.compile(r'Br[+-]?|Cl[+-]?|H[0-9]?[+-]?|[a-zA-Z][+-]?|[=#]|[0-9]{1}|[()]')
AMINO_ACID_REGEX = re.compile(r'[nN]H[23]?\+')
