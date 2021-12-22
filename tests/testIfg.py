""" Use this file to test a single SMILES code for its functional groups and check the output of values

"""
from helpers import createFgDataDict
from ifg import ifg
import pandas

# User Test Molecule
SMILES = "OC(=O)CC1CCC(=O)CC1"
REFCODE = "KUZQIG"
###########################

def print_data(matches_list):
    """ Displays data in allFgs or preciseFgs via pandas
        matches_list (list):  either allFgs or preciseFgs
    """
    # Organize and format data for pandas dataframe inputs
    data = createFgDataDict(matches_list)
    fg_name_to_template = {f.NAME: f.SMILES for f in matches_list }
    columns = ["Name", "Template", "Count"]
    dat = [[k, fg_name_to_template[k], v] for (k,v) in data.items()]

    # Create and display pandas dataframe
    df = pandas.DataFrame(columns=columns, data=dat)
    print(df)


# Run IFG and collect matches_list output
fgs = ifg(SMILES, REFCODE)
print("\nOUTPUT \n")
print(fgs)
print("\nAll fgs")
print_data(fgs.allFgs)
print("\nPrecise fgs")
print_data(fgs.preciseFgs)
print("\n")