"""Configuration file to new SMILES codes lists, to the FGlist, and to the generated excel output file
"""

from pathlib import Path

# Directory paths
root_path = Path('.')
src_path = Path("src")
resources_path = Path("resources")

# Test paths and view files in directories
# [x for x in root_path.iterdir()]
# [x for x in src_path.iterdir()]
# [x for x in resources_path.iterdir()]

# SMILES Codes path
SMILESPATH = resources_path / 'smiles.txt'              # Change target SMILES file
# for x in open(SMILESPATH.resolve(), 'r'):             # Test SMILES path
#     print(x)

# Functional Groups path
FGSPATH = resources_path / 'FGlist.txt'             # Change target FG list file
# for x in open(FGSPATH.resolve(), 'r'):              # Test FG path
#     print(x)

# Output file path
out_file = 'output/FunctionalGroups.xlsx'

