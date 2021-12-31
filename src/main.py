""" Functional script that loops over the target set of smiles codes and determines thier functional groups 
"""

from ifg import ifg
import pandas as pd
import re
from helpers import createFgDataDict
import numpy as np
from collections import defaultdict
from config import FGSPATH, SMILESPATH
from tqdm import tqdm

def main():
    """ Returns a dictionary of two dataframes containing the data about functional groups for a given set of smiles codes

        Notes:
            allDf is the dataframe where "all" functional groups are counted
            preciseDf is the dataframe where the "precise" functional groups are counted

    """


    FGlist_data = [                                           # Fetch the (smiles, name) functional group pairs from FGlist.txt
        [y.strip() for y in x.split(' ')]                     # Split into a formatted (smiles, name) pair
        for x                                                 # Done for each line
        in open(FGSPATH.resolve(), "r+").readlines()          # Create list using all lines in FGlist.txt
    ]
    FGlist_data.append(["OH", "Alcohol"])
    FGlist_data.append(["NR", "PrimaryAmine"])

    FGnames = sorted({name[1] for name in FGlist_data})        # FGnames in a sorted set

    cyclicGroups = ["Cyclic" + group for group in FGnames]      # Cyclic nomenclatures for FGnames
    aromaticGroups = ["Aromatic" + group for group in FGnames]  # Aromatic nomenclatures for FGnames

    columns = [                                                 # Combine calssified and base functional group names into a single list of strings (no tuples)
        name for names                                          # Loop over tuples created by zip
        in zip(FGnames, cyclicGroups, aromaticGroups)           # Tuple of an FG with its cyclic and aromatic nomenclature, same bases
        for name in names                                       # Loop over elements within tuple 
    ]

    columns.insert(0, "SMILES")
    columns.insert(0, "Refcode")

    properties = ['aromaticRingCount',                          # Extra data ontop of base functional group analysis
                  'nonAromaticRingCount',
                  'ringCount',
                  'totalAlcohols',
                  'AminoAcid',
                  'Bromine',
                  'Chlorine',
                  'Iodine',
                  'Fluorine'
    ]

    for prop in properties:                                     # Attach additional prop names to column list
        columns.append(prop)
    
    # Set up data lists
    allDfData = []
    preciseDfData = []

    # Set up performance metric varbiales for execution times
    p = 0
    z = 0

    with tqdm(total=len(open(SMILESPATH.resolve()).readlines())) as bar:
        for text in open(SMILESPATH.resolve()):                     # Loop over all smiles codes in this current dataset

            line = re.compile(r'\S+').findall(text)                 # Get line from smiles text list
            (temp, smiles, refcode) = line                                # Extract the line into variables
            functionalGroups = ifg(smiles, refcode)                 # Determine the functional groups based on the input SMILES code

            p+=functionalGroups.dfs_time.total_seconds()
            z+=functionalGroups.filter_time.total_seconds()

            propData = {                                            # Add additional properties baesd on Molecule
                **functionalGroups.ringData,
                "totalAlcohols": len(functionalGroups.ALCOHOLICINDICES),
                "AminoAcid": "Yes" if functionalGroups.AMINOACID else "No",
                **functionalGroups.HALOGENS
            }

            all_fgs = defaultdict(int, {                             # Combine properties with all_fgs to get all_fgs dict
                **createFgDataDict(functionalGroups.all_fgs),    
                **propData
            })

            exact_fgs = defaultdict(int, {                         # Combine properties with exact_fgs to get exact_fgs dict
                **createFgDataDict(functionalGroups.exact_fgs),
                **propData
            })
                                                                    # Decompse dictionary into values only list which mirrors the 
                                                                    # index positioning of the columns in the all dataframe
            ALL_DATA = [                                             
                all_fgs[name] if all_fgs[name]
                else np.nan
                for name in columns[2:]
            ]
                                                                    # Do this for precise FGs as well
            EXACT_DATA = [
                exact_fgs[name] if exact_fgs[name]
                else np.nan
                for name in columns[2:]
            ]

            for _id in [smiles, refcode]:                           # Prepend smiles and refcode to satify structure of columns
                ALL_DATA.insert(0, _id)
                EXACT_DATA.insert(0, _id)

            allDfData.append(ALL_DATA)                            # Locate a new row indexed by refcode
            preciseDfData.append(EXACT_DATA)                        # For both dataframes

            bar.update(1)                                           # Increment the progress bar once the Molecule is processed


    allDf = pd.DataFrame(columns=columns, data=allDfData)                       # Initialize the all data frame with column names
    preciseDf = pd.DataFrame(columns=columns, data=preciseDfData)                   # And the precise data frame with column names

    allDf.set_index("Refcode")                                  # Index the DF's by REFCODE
    preciseDf.set_index("Refcode")

    dfs = {}                                                    # Return nan-filtered dataframes (i.e. columns with all NaN are removed)
    dfs.update({"allDf": allDf.dropna(axis=1, how='all')})
    dfs.update({"preciseDf": preciseDf.dropna(axis=1, how='all')})


    # Show the total performance
    total = p + z
    print("DFS Performance Evaluation")
    print("all seconds = ", p)
    print("percentage dfs = ", (p/total)*100)
    print("average time = ", p/831) # Change to varaible length
    print("\n")

    print("FILTERING Performance Evaluation")
    print("all seconds = ", z)
    print("percentage filtering = ",(z/total)*100)
    print("average time = ", z/831) # Change to varaible length
    print("\n")


    return dfs
