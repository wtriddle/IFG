from pathlib import Path

root_path = Path('.')
src_path = Path("src")
resources_path = Path("resources")

# Test paths and view files in directories
[x for x in root_path.iterdir()]
[x for x in src_path.iterdir()]
[x for x in resources_path.iterdir()]

SMILESPATH = resources_path / 'smiles.txt'          # Change target SMILES file
SMILESPATH.resolve()
for x in open(SMILESPATH.resolve(), 'r'):
    print(x)

FGSPATH = resources_path / 'FGlist.txt'             # Change target FG list file
for x in open(FGSPATH.resolve(), 'r'):
    print(x)



