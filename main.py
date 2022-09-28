"""Creates functional group excel data sheet given a target smiles code set from a text file"""

import pandas
from tqdm import tqdm
from molecule import Molecule
from datetime import datetime
from config import STRUCTURES_PATH, MAIN_OUTPUT_PATH

def main():
    """Creates a chemical data excel sheet from a set of SMILES according to the data extracted from the Molecule class"""
    
    ##### Execution Time Variable #####
    now = datetime.now()

    ##### Data Variables #####
    all_data: list[dict] = []
    exact_data: list[dict] = []

    ##### Input Structure Data Load #####
    STRUCTURES = open(str(STRUCTURES_PATH.resolve()), "r+").readlines()

    ##### Structure Bar Status #####
    with tqdm(total=len(STRUCTURES)) as bar:

        ##### SMILES Structure Loop #####
        for (_, smiles, refcode) in [[y for y in x.strip().split(' ') if y] for x in STRUCTURES]:

            ##### Molecule Data #####
            mol = Molecule(smiles, refcode, type='mol')

            ##### All Functional Group Format Data #####
            all_data.append({
                "Refcode": mol.name,
                "SMILES": smiles,
                "Aromatic Rings": mol.aromatic_ring_count,
                "Non Aromatic Rings": mol.non_aromatic_ring_count,
                "Rings": mol.total_ring_count,
                "AminoAcid": "Yes" if mol.amino_acid else "No",
                **mol.functional_groups_all,
            })

            ##### Exact Functional Group Format Data #####
            exact_data.append({
                "Refcode": mol.name,
                "SMILES": smiles,
                "Aromatic Rings": mol.aromatic_ring_count,
                "Non Aromatic Rings": mol.non_aromatic_ring_count,
                "Rings": mol.total_ring_count,
                "AminoAcid": "Yes" if mol.amino_acid else "No",
                **mol.functional_groups_exact,
            })

            ##### Status Bar Update #####
            bar.update(1)                                               # Increment the progress bar once smiles finishes processing

    ##### Execution Time Evaluation #####
    print("execution time = ", datetime.now() - now, "s")
    print("\n")

    ##### Pandas Dataframe #####
    df_all = pandas.DataFrame(all_data).fillna(0).set_index("Refcode")
    df_exact = pandas.DataFrame(exact_data).fillna(0).set_index("Refcode")

    ##### Excel Exporter #####
    writer = pandas.ExcelWriter(str(MAIN_OUTPUT_PATH.resolve()), engine="xlsxwriter")

    ##### All Functional Groups Data Sheet Export #####
    df_all.to_excel(writer, sheet_name="all_data", freeze_panes=(1, 1))
    all_sheet = writer.sheets["all_data"]
    all_sheet.set_column(0, 0, 13)
    for i in range(0, len(df_all.columns)):
        all_sheet.set_column(i+1, i+1, len(df_all.columns[i])+7)

    ##### Exact Functional Groups Data Sheet Export #####
    df_exact.to_excel(writer, sheet_name="exact_data", freeze_panes=(1, 1))
    exact_sheet = writer.sheets["exact_data"]
    exact_sheet.set_column(0, 0, 13)
    for i in range(0, len(df_exact.columns)):
        exact_sheet.set_column(i+1, i+1, len(df_exact.columns[i])+7)

    ##### Excel File Save #####
    writer.save()

if __name__ == "__main__":
    main()