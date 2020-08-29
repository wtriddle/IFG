import pandas as pd
from ifgTest import ifgTest
import re
from collections import Counter
from helpers import createFgDataDict
import sys
import numpy as np
from collections import defaultdict
import os


def identifyFunctionalGroups(allSheet, preciseSheet, verboseMode):
    """
        Algorithm to create the dataframe object that contains refcode to functional group relation
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
    properties = ['aromaticRings',
                  'nonAromaticRings',
                  'ringCount',
                  'totalAlcohols',
                  'AminoAcid']

    for prop in properties:
        columns.append(prop)

    # Dataframe containers for SMILES (index) to Functioanl group (columns) counts relation
    allDf = pd.DataFrame(columns=columns)
    preciseDf = pd.DataFrame(columns=columns)

    # Loop over all smiles codes in this current dataset
    for i, text in enumerate(open(os.getcwd() + '/src/resources/smiles.txt')):

        # Retrive Refcode and its smiles representation from set of smiles codes
        line = re.compile(r'\S+').findall(text)
        (smilesAlt, smiles, refcode) = line

        # Process the smiles code with the ifgTest algorithm and retireve output data in the form of a dictionary
        # Inherets molecule properties from Molecule class
        functionalGroups = ifgTest(smiles, refcode)

        if verboseMode:
            print(refcode, " : ", smiles)

        # Retrieve molecule data from functionalgroups algorithm and molecule data
        propData = {
            **functionalGroups.RINGDICT,
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
            print(columns)

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

        # Inser smiles and refcode to satify structure of columns
        for _id in [smiles, refcode]:
            allData.insert(0, _id)
            preciseData.insert(0, _id)

        # Locate a new row indexed by refcode
        if allSheet:
            allDf.loc[refcode] = allData

        if preciseSheet:
            preciseDf.loc[refcode] = preciseData

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
