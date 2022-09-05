"""File-System Path Configurations"""

from pathlib import Path

##### Target Molecular SMILES Codes #####
STRUCTURES_PATH = Path("resources") / 'smiles.txt'

##### Target Functional Groups #####
FUNCTIONAL_GROUPS_PATH = Path("resources") / 'FGlist.txt'

##### Target Main Output Excel Sheet #####
MAIN_OUTPUT_PATH = Path("output") / 'output.xlsx'

##### Target Analysis Output Excel Sheet #####
ANALYSIS_OUTPUT_PATH = Path("output") / 'stats.xlsx'

##### Target Structure Bandgaps Input Excel Sheet #####
BANDGAPS_PATH = Path("resources") / 'CrystalData.xlsx'