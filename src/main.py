""" Module containing the functional script that loops over a CHON set of smiles codes and determines thier functional groups """

from ifg import ifg
import pandas as pd
import re
from helpers import createFgDataDict
import sys
import numpy as np
from collections import defaultdict
import os
from progress.bar import IncrementalBar


def identifyFunctionalGroups(allSheet, preciseSheet, verboseMode):
    """ Returns a dictionary of two dataframes containing the data about functional groups for the entire smiles code set 


        allSheet (bool) : Determines if script is to return functional group data in the all format
        preciseSheet (bool) : Determines if script is to return functional group data in the precise format    
        verboseMode (bool) : Determines to print functioanl group contents out to string or keep the same

        Notes:
            Dataframes are allDf and preciseDf, corellating to what type of functional group data it contains.

    """

    # fgNames is a list of strings with each functional group name that participated in this analysis
    # Retrieve functional group names (strings) from FGlist text file
    fgNames = pd.read_csv(os.getcwd() + '/src/resources/FGlist.txt',
                          sep=" ", header=None)
    fgNames.columns = ['template', 'name']

    # Alcohol/Primary amine manuel entry b/c no template in FGlist.txt
    fgNames = fgNames.append(
        pd.Series(["OH", "Alcohol"], index=fgNames.columns),
        ignore_index=True
    )
    fgNames = fgNames.append(
        pd.Series(["NR", "PrimaryAmine"], index=fgNames.columns),
        ignore_index=True
    )
    fgNames = sorted({name for name in fgNames['name']})

    # Additional cyclic/aromaic nomenclatures for all functional groups in base string list
    cyclicGroups = ["Cyclic" + group for group in fgNames]
    aromaticGroups = ["Aromatic" + group for group in fgNames]

    # Combine calssified and base functional group names into a single list
    # (each item is of type string)
    columns = [
        name for names
        in zip(fgNames, cyclicGroups, aromaticGroups)
        for name in names
    ]

    # Insert extra string info
    columns.insert(0, "SMILES")
    columns.insert(0, "Refcode")

    # Extra info ontop of base functional group analysis
    properties = ['aromaticRingCount',
                  'nonAromaticRingCount',
                  'ringCount',
                  'totalAlcohols',
                  'AminoAcid']

    for prop in properties:
        columns.append(prop)

    if verboseMode:
        print("Created column names")

    # Dataframe containers for SMILES (index) to Functioanl group (columns) counts relation
    allDf = pd.DataFrame(columns=columns)
    preciseDf = pd.DataFrame(columns=columns)

    # Progress bar for script execution, goes to number of smiles codes
    bar = IncrementalBar('Anlalyzing smiles codes', max=len(
        open(os.getcwd() + '/src/resources/smiles.txt').readlines()
    ))

    # Loop over all smiles codes in this current dataset
    for i, text in enumerate(open(os.getcwd() + '/src/resources/smiles.txt')):

        # Retrive Refcode and its smiles representation from set of smiles codes
        line = re.compile(r'\S+').findall(text)
        (smilesAlt, smiles, refcode) = line

        # Process the smiles code with the ifg algorithm and retireve output data in the form of a dictionary
        # Inherets molecule properties from Molecule class
        functionalGroups = ifg(smiles, refcode)

        if verboseMode:
            print(refcode, " : ", smiles)

        # Retrieve molecule data from functionalgroups algorithm and molecule data
        propData = {
            **functionalGroups.ringData,
            "totalAlcohols": len(functionalGroups.ALCOHOLICINDICES),
            "AminoAcid": "Yes" if functionalGroups.AMINOACID else "No"
        }

        allFgs = defaultdict(int, {
            **createFgDataDict(functionalGroups.allFgs),
            **propData
        })

        preciseFgs = defaultdict(int, {
            **createFgDataDict(functionalGroups.preciseFgs),
            **propData
        })

        # Print results of each algorithm if verbose mode is turned on
        if verboseMode:
            print(allFgs)
            print(preciseFgs)

        # Loop over all possible FGs and,
        # in the order of columns as they appear in dataframe,
        # insert the count of a certain FG to the named column slot,
        # if the name was not in the dictionary, NaN of that FG were found
        allData = [
            allFgs[name] if allFgs[name]
            else np.nan
            for name in columns[2:]
        ]

        preciseData = [
            preciseFgs[name] if preciseFgs[name]
            else np.nan
            for name in columns[2:]
        ]

        # Insert smiles and refcode to satify structure of columns
        for _id in [smiles, refcode]:
            allData.insert(0, _id)
            preciseData.insert(0, _id)

        # Locate a new row indexed by refcode
        if allSheet:
            allDf.loc[refcode] = allData

        if preciseSheet:
            preciseDf.loc[refcode] = preciseData

        bar.next()

    bar.finish()

    # Set index by Refcode column
    allDf.set_index("Refcode")
    preciseDf.set_index("Refcode")

    # Return nan-filtered dataframes
    # If a column is entirley Nan, it is dropped
    dfs = {}
    if allSheet:
        dfs.update({"allDf": allDf.dropna(axis=1, how='all')})
    if preciseSheet:
        dfs.update({"preciseDf": preciseDf.dropna(axis=1, how='all')})

    return dfs
