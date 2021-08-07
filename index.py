""" Top level wrapper script that runs main() from the src folder. 

    To run this script, open a terminal and run:
    python3 index.py

    Notes:
        To run this script, the PYTHONPATH variable must be set to include the /src file as a target:

        On Linux:
        1. Go to ~/.bashrc
        2. add the following line:
        export PYTHONPATH="(/PATH/TO/IFG/GOES/HERE)$PATH"
        example:
        export PYTHONPATH="/media/wtrid/Wills_2nd_drive/Theses/CHMSeniorThesis/Resources/Code/IFG:$PATH"

        On Windows:
        1. Go to environment variables
        2. Click on add new environment variables
        3. Create a new environemnt variable called PYTHONPATH
        4. Add its value as the path to the top level IFG folder
        example:
        C:/Users/wtriddle/OneDrive/Desktop/IFG

"""

from main import main 
import pandas as pd
import json
from config import out_file


def index():
    """ Top level wrapper script to handle excel sheet outputs and formatting of columns """

    data = main()                                            # Retrieve the data from IFG
    writer = pd.ExcelWriter(out_file, engine="xlsxwriter")   # Create pandas excel writer with xlsxwriter as engine

    try:
        data['allDf'].to_excel(
            writer,
            sheet_name="All Functional Groups",
            na_rep=0,                                       
            freeze_panes=(1, 1)                             # Freeze columns and rows
        )
        allSheet = writer.sheets['All Functional Groups']   # XlsxWriter sheet object
        allSheet.set_column(1,                              # Set sheet column width based on longest name
                            len(data['allDf'].columns), 
                            21, 
                            None
        )

        allJsonData = data['allDf'].to_json(orient="index")     # JSON format of output data
        d = json.loads(allJsonData)
        with open("allData.json", "w", encoding="utf-8") as f:  # Save JSON data as allData.json
            json.dump(d, f, indent=4)                           
    except:
        pass

    try:
        data['preciseDf'].to_excel(
            writer,
            sheet_name="Precise Functional Groups",
            na_rep=0,
            freeze_panes=(1, 1)                                     # Freeze columns and rows
        )
        preciseSheet = writer.sheets['Precise Functional Groups']   # XlsxWriter sheet object
        preciseSheet.set_column(1,                                  # Set sheet column width based on longest name
                                len(data['preciseDf'].columns), 
                                21, 
                                None
        )

        preciseJsonData = data['preciseDf'].to_json(orient="index") # JSON format of data
        d = json.loads(preciseJsonData)
        with open("preciseData.json", "w", encoding="utf-8") as f:  # Save JSON data as preciseData.json
            json.dump(d, f, indent=4)
    except:
        pass

    try:
        writer.save()                                               # Save the sheet
    except:
        raise EnvironmentError("File could not be saved")


if __name__ == '__main__':
    index()
