""" Use this file to test a single SMILES code for its functional groups and check the output of values

"""
from helpers import createFgDataDict
from ifg import ifg
import pandas

# User Test Molecule
SMILES = "COC(=O)C(C)(O)NC(C)=O"
REFCODE = "LAQSOM"
###########################

def print_data(matches_list):
    """ Displays data in all_fgs or exact_fgs via pandas
        matches_list (list):  either all_fgs or exact_fgs
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
print_data(fgs.all_fgs)
print("\nExact fgs")
print_data(fgs.exact_fgs)
print("\n")