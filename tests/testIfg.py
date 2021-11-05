""" Use this file to test a single SMILES code for its functional groups and check the output of values

"""

from re import template
from helpers import createFgDataDict
from ifg import ifg
import pandas

fgs = ifg("O=N(=O)c1ccc(NN=C2C3CC(C=C3)C3CCCC23)c(c1)N(=O)=O",  "APYFEB01")
data = createFgDataDict(fgs.preciseFgs)
print(data)
dict_to_template = {}
for f in fgs.preciseFgs:
    dict_to_template[f.NAME] = f.SMILES
print("Name\t Template\t Count")
dat = []
for  i, (k, v) in enumerate(data.items()):
    template = dict_to_template[k]
    dat.append([])
    dat[i] = [k, template, v]
dat
df = pandas.DataFrame(columns=["Name", "Template", "Count"], data=dat)
df

for d in data:
    print(fgs.preciseFgs)
for f in fgs.preciseFgs:
    print(f.SMILES)
print("\n".join("{}\t{}".format(k, v) for k, v in data.items()))