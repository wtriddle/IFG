""" Handles foramtting of inputs and outputs of IFG to retireve the data and export it to two dataframes
    Gives performance report of computation time 
    Notes:
        **dict notation inside of dictionary construction dumps 
        all key-value pairs from **dict into new dict

        Ex:
        newDict = {
            **otherDict,    (All k-v pairs dumped into newDict)
            key: value      (New k-v pair into newDict as well)
        }
"""

from ifg import ifg
import pandas as pd
from helpers import createFgDataDict
import numpy as np
from collections import defaultdict
from config import FGSPATH, SMILESPATH
from tqdm import tqdm

def main():
    """ Returns a dictionary of two dataframes containing the data about functional groups for a given set of smiles codes
        Gives performance report of IFG (DFS process and Filtering processes included)
        Notes:
            allDf is the dataframe where "all" functional groups are counted
            precise is the dataframe which does not use the overlapping filter
    """

    ##### Functional Group SMILES & NAME Loading #####
    FGlist_data = [                                             # Fetch the (smiles, name) functional group pairs from FGlist.txt
        [y.strip() for y in x.split(' ')]                       # Split into a formatted (smiles, name) pair
        for x                                                   # Done for each line
        in open(FGSPATH.resolve(), "r+").readlines()            # Create list using all lines in FGlist.txt
    ]
    FGlist_data.append(["OH", "Alcohol"])                       # Alochol not identified using FGlist.txt
    FGlist_data.append(["NR", "PrimaryAmine"])                  # PrimaryAmine not identified using FGlist.txt

    ##### Functional Group Dataframe Column NAME Creation #####
    FGnames = sorted({name[1] for name in FGlist_data})         # Names of all identifiable functional groups by IFG
    cyclicGroups = ["Cyclic" + group for group in FGnames]      # Cyclic nomenclatures for FGnames
    aromaticGroups = ["Aromatic" + group for group in FGnames]  # Aromatic nomenclatures for FGnames

    columns = [                                                 # List of all functional group named columns in exported Dataframe from IFG
        name for names                                          # Loop over tuples created by zip
        in zip(FGnames, cyclicGroups, aromaticGroups)           # Tuple list in format like: [("Ether", "CyclicEther", "AromaticEther"), ...]
        for name in names                                       # Loop over elements within tuple list to obtain all FG names with aromatic, non-aromatic, non-cyclic nomenclature
    ]

    ##### Structure NAMES & Molecular Data NAMES Column Creations #####
    columns.insert(0, "SMILES")
    columns.insert(0, "Refcode")
    for prop in [ 'aromaticRingCount',           # Molecular based SMILES data (properties) that comes with IFG (from SSCD) on SMILES code
                  'nonAromaticRingCount',
                  'ringCount',
                  'totalAlcohols',
                  'AminoAcid',
                  'Bromine',
                  'Chlorine',
                  'Iodine',
                  'Fluorine'
    ]:                                     # Add each Molecular data NAME to column names
        columns.append(prop)
    
    ##### I/O and Performance Metric Variables #####
    allDfData = []
    preciseDfData = []

    dfs_time_total = 0
    filter_time_total = 0

    SMILES_FILE = open(SMILESPATH.resolve(), "r+")
    SMILES_LINES = SMILES_FILE.readlines()
    TOTAL_STRUCTURES_NUM = len(SMILES_LINES)

    ##### SMILES Set Loop & IFG Data Collection #####
    with tqdm(total=TOTAL_STRUCTURES_NUM) as bar:           # Show progress of smiles processing

        for (temp, smiles, refcode) in [                    # Tuple split SMILES string and structure refcode (adjusted to fit format of target file)
            [y.strip() for y in x.split(' ') if y != '']               # Split each space-divided string into a formatted tuple set
            for x                                           # Do this for each line
            in SMILES_LINES                                 # Loop is over all SMILES strings in selected SMILES list for data collection 
        ]:

            functionalGroups = ifg(smiles, refcode)                 # Collect molecular and functional group data from the SMILES code

            dfs_time_total+=functionalGroups.dfs_time.total_seconds()
            filter_time_total+=functionalGroups.filter_time.total_seconds()

            molData = {                                             # Retrieve molecular data (SSCD from smiles)
                **functionalGroups.ringData,
                "totalAlcohols": len(functionalGroups.ALCOHOLICINDICES),
                "AminoAcid": "Yes" if functionalGroups.AMINOACID else "No",
                **functionalGroups.HALOGENS
            }

            all_fgs = defaultdict(int, {                           # Combine all_fgs and molecular data for all_fgs format
                **createFgDataDict(functionalGroups.all_fgs),    
                **molData
            })

            exact_fgs = defaultdict(int, {                         # Combine exact_fgs and molecular data for exact_fgs format
                **createFgDataDict(functionalGroups.exact_fgs),
                **molData
            })

            ALL_DATA = [                                            # Create row of ALL_DATA dataframe based on column positioning
                all_fgs[name] if all_fgs[name]                      # Dictionary key from name
                else np.nan
                for name in columns[2:]                             # Name in list order of columns so it matches row in dataframe with column positions
            ]
                                                                    
            EXACT_DATA = [                                          # Create row of ALL_DATA dataframe based on column positioning
                exact_fgs[name] if exact_fgs[name]
                else np.nan
                for name in columns[2:]
            ]

            for _id in [smiles, refcode]:                           # Prepend smiles and refcode to identify structure
                ALL_DATA.insert(0, _id)
                EXACT_DATA.insert(0, _id)

            allDfData.append(ALL_DATA)                              # Append the ALL_DATA row allDf data
            preciseDfData.append(EXACT_DATA)                        # Append the EXACT_DATA row preciseDf data

            bar.update(1)                                           # Increment the progress bar once smiles finishes processing

    ##### IFG Performance Report #####
    total_IFG = dfs_time_total + filter_time_total
    print("DFS Performance Evaluation")
    print("all seconds = ", dfs_time_total)
    print("percentage dfs = ", (dfs_time_total/total_IFG)*100)
    print("average time = ", dfs_time_total/TOTAL_STRUCTURES_NUM)                # Change to varaible length
    print("\n")

    print("FILTERING Performance Evaluation")
    print("all seconds = ", filter_time_total)
    print("percentage filtering = ",(filter_time_total/total_IFG)*100)
    print("average time = ", filter_time_total/TOTAL_STRUCTURES_NUM)             # Change to varaible length
    print("\n")

    ##### Data Frame Creation and Export #####
    allDf = pd.DataFrame(columns=columns, data=allDfData)                 
    preciseDf = pd.DataFrame(columns=columns, data=preciseDfData)        

    allDf.set_index("Refcode")
    preciseDf.set_index("Refcode")

    dfs = {}
    dfs.update({"allDf": allDf.dropna(axis=1, how='all')})          # Add the allDf to exported data from main
    dfs.update({"preciseDf": preciseDf.dropna(axis=1, how='all')})  # Add the preciseDf to exported data from main

    return dfs
