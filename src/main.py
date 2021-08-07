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

    FGtemps = pd.read_csv(FGSPATH.resolve(), sep=" ", header=None)  # Dataframe of functional group templates and names
    FGtemps.columns = ['template', 'name']                          # Columns for FGlist

                                                                # No Alcohol or Primary Amine in FGlist, manually determined
    FGtemps = FGtemps.append(                                   # Manual Alcohol template entry
        pd.Series(
            ["OH", "Alcohol"], index=FGtemps.columns
        ),
        ignore_index=True
    )
    FGtemps = FGtemps.append(                                   # Manual Primary Amine template entry
        pd.Series(
            ["NR", "PrimaryAmine"], index=FGtemps.columns
        ),
        ignore_index=True
    )
    FGnames = sorted({name for name in FGtemps['name']})        # FGnames in a sorted set

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

    allDf = pd.DataFrame(columns=columns)                       # Initialize the all data frame with column names
    preciseDf = pd.DataFrame(columns=columns)                   # And the precise data frame with column names

    with tqdm(total=len(open(SMILESPATH.resolve()).readlines())) as bar:
        for text in open(SMILESPATH.resolve()):                     # Loop over all smiles codes in this current dataset

            line = re.compile(r'\S+').findall(text)                 # Get line from smiles text list
            (smiles, refcode) = line                                # Extract the line into variables
            functionalGroups = ifg(smiles, refcode)                 # Determine the functional groups based on the input SMILES code

            propData = {                                            # Add additional properties baesd on Molecule
                **functionalGroups.ringData,
                "totalAlcohols": len(functionalGroups.ALCOHOLICINDICES),
                "AminoAcid": "Yes" if functionalGroups.AMINOACID else "No",
                **functionalGroups.HALOGENS
            }

            allFgs = defaultdict(int, {                             # Combine properties with allFgs to get allFGs dict
                **createFgDataDict(functionalGroups.allFgs),    
                **propData
            })

            preciseFgs = defaultdict(int, {                         # Combine properties with preciseFgs to get preciseFGs dict
                **createFgDataDict(functionalGroups.preciseFgs),
                **propData
            })
                                                                    # Decompse dictionary into values only list which mirrors the 
                                                                    # index positioning of the columns in the all dataframe
            allData = [                                             
                allFgs[name] if allFgs[name]
                else np.nan
                for name in columns[2:]
            ]
                                                                    # Do this for precise FGs as well
            preciseData = [
                preciseFgs[name] if preciseFgs[name]
                else np.nan
                for name in columns[2:]
            ]

            for _id in [smiles, refcode]:                           # Prepend smiles and refcode to satify structure of columns
                allData.insert(0, _id)
                preciseData.insert(0, _id)

            allDf.loc[refcode] = allData                            # Locate a new row indexed by refcode
            preciseDf.loc[refcode] = preciseData                    # For both dataframes

            bar.update(1)                                           # Increment the progress bar once the Molecule is processed

    allDf.set_index("Refcode")                                  # Index the DF's by REFCODE
    preciseDf.set_index("Refcode")

    dfs = {}                                                    # Return nan-filtered dataframes (i.e. columns with all NaN are removed)
    dfs.update({"allDf": allDf.dropna(axis=1, how='all')})
    dfs.update({"preciseDf": preciseDf.dropna(axis=1, how='all')})

    return dfs
