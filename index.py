""" Top level script for running IFG and exporting excel sheet data of structure-functional group data

    SETUP & EXECUTION:

        - The PYTHONPATH variable must be set to include the /src file as a target to find required files for execution:

        On Linux:
        1. Go to ~/.bashrc
        2. add the following line:
        export PYTHONPATH="(/PATH/TO/IFG/GOES/HERE)$PATH"
        example:
        export PYTHONPATH="/media/wtrid/Wills_2nd_drive/Theses/CHMSeniorThesis/Resources/Code/IFG:$PATH"
        3. cd ~/IFG
        4. python index.py
            ...
        5. View results in targeted output folder 

        On Windows:
        1. Go to environment variables
        2. Click on add new environment variables
        3. Create a new environemnt variable called PYTHONPATH
        4. Add its value as the path to the top level IFG folder
            example:
            C:/Users/wtriddle/OneDrive/Desktop/IFG
        5. Restart computer
        6. Run this script in the anaconda base environment with:
            C:/Users/wtriddle/OneDrive/Desktop/IFG> conda activate base
            (base) C:/Users/wtriddle/OneDrive/Desktop/IFG> python ./index.py
            ...
        7. View results in targeted output folder

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
