""" Top level script run with py RunIFG.py [args, [...]] """

""" Create a dependencies list so people can install environemnt right off the bat
    If the script has run once fully and failed at the end, add functionality to load the .json files directly to excel instead
    Of re-calculating the results
"""

import sys
import getopt
from main import identifyFunctionalGroups
import pandas as pd
import pathlib
import os
import json


def main(argv):
    preciseSheet = False
    allSheet = False
    verbose = False
    file = 'output/FunctionalGroups.xlsx'
    try:
        args, opts = getopt.getopt(
            argv, 'aphfv', ["all", "precise", "help", "verbose", "file="])
    except getopt.GetoptError:
        print('Bad usage. Type smiles.py -h for help')
        sys.exit()
    for opt, arg in args:
        if opt in ['-h', '--help']:
            # Help option
            print(
                """smiles.py [ -a | --all <bool> ] [ -p | --precise <bool> ] [ -v || --verbose <bool> ] [-h | --help] [ -f | --file <path>]
                """)
            sys.exit()
        elif opt in ['-a', '--all']:
            allSheet = True
        elif opt in ['-p', '--precise']:
            preciseSheet = True
        elif opt in ['-v', '--verbose']:
            verbose = True
        elif opt in ['-f', '--file']:
            file = 'output/' + arg + '.xlsx'
    if not args:
        preciseSheet = True
        allSheet = True
        file = 'output/FunctionalGroups.xlsx'

    # Retrieve functional groups data, returns a dictionary
    data = identifyFunctionalGroups(allSheet, preciseSheet, verbose)

    # Create pandas excel writer with xlsxwriter as engine
    writer = pd.ExcelWriter(file, engine="xlsxwriter")

    # Input respective data into excel sheets if they were created
    try:
        data['allDf'].to_excel(
            writer,
            sheet_name="All Functional Groups",
            na_rep=0,
            freeze_panes=(1, 1)
        )
        # XlsxWriter sheet object
        allSheet = writer.sheets['All Functional Groups']
        # Set sheet column width
        allSheet.set_column(1, len(data['allDf'].columns), 21, None)

        # JSON format of data
        allJsonData = data['allDf'].to_json(orient="index")
        d = json.loads(allJsonData)
        with open("allData.json", "w", encoding="utf-8") as f:
            json.dump(d, f, indent=4)
    except:
        pass

    try:
        data['preciseDf'].to_excel(
            writer,
            sheet_name="Precise Functional Groups",
            na_rep=0,
            freeze_panes=(1, 1)
        )
        # XlsxWriter sheet object
        preciseSheet = writer.sheets['Precise Functional Groups']
        # Set sheet column width
        preciseSheet.set_column(1, len(data['preciseDf'].columns), 21, None)

        # JSON format of data
        preciseJsonData = data['preciseDf'].to_json(orient="index")
        d = json.loads(preciseJsonData)
        with open("preciseData.json", "w", encoding="utf-8") as f:
            json.dump(d, f, indent=4)
    except:
        pass

    # Save the sheet
    try:
        writer.save()
    except:
        raise EnvironmentError("File could not be saved")


if __name__ == '__main__':
    main(sys.argv[1:])
