"""Analyzes a functional group dataset and a bandgap dataset for quantitative structure-property relationships"""

from datetime import datetime
from itertools import chain, combinations
import math
import os

import numpy
import pandas
from progress.spinner import PieSpinner

##### Spinner Progress #####
spinner = PieSpinner("Processing ")
spinner.next()

##### Target Main Output Excel Sheet #####
MAIN_OUTPUT_PATH = os.path.dirname(__file__) + '/output/functional_groups.xlsx'
"""Excel file generated by the main.py script"""

##### Target Analysis Output Excel Sheet #####
ANALYSIS_OUTPUT_PATH = os.path.dirname(__file__) + '/output/stats.xlsx'
"""Excel file generated by this script"""

##### Target Structure Bandgaps Input Excel Sheet #####
BANDGAPS_PATH = os.path.dirname(__file__) + '/output/CrystalData.xlsx'
"""Excel file with the bandgaps used by this script"""

##### User Defined combinational Functional Group Sets #####
SETS = [
    ["PrimaryAmine", "SecondaryAmine", "TertiaryAmine"],
]

##### combinational Functional Groups Color Schemes #####
class Color():
    def __init__(self, label_color: str, bar_color: str):
        self.label_color: str = label_color
        self.bar_color: str = bar_color

COLORS = [
    Color(label_color, bar_color)
    for (label_color, bar_color) in
    [
        ("#ffffff", "#0b84a5"),
        ("#000000", "#f6c85f"),
        ("#ffffff", "#6f4e7c"),
        ("#000000", "#9dd866"),
        ("#ffffff", "#ca472f"),
        ("#000000", "#ffa056"),
        ("#000000", "#8dddd0"),
    ]
]

##### Excel Sheet Target Output File With Writer Engine #####
writer = pandas.ExcelWriter(ANALYSIS_OUTPUT_PATH)   # type: ignore

##### Powerset Function #####
def powerset(iterable):
    """Produce all subsets, including the full set, of a given iterable"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


##### Combinational Query Generator Function #####
def generateCombinations(fg_set: "list[str]") -> "list[str]":
    """Generate all combinational functional group queries for a pandas dataframe using the functional group names listed in fg_set"""
    
    return [
        " & ".join(
            [
                "`" + fg + "` > 0" if fg in non_zero_fgs 
                else "`" + fg + "` == 0" 
                for fg in fg_set
            ]
        )
        for non_zero_fgs in list(powerset(fg_set)) if non_zero_fgs
    ]


##### Column-Based Functional Group Statistics Generator Function #####
def createStatisticalMetrics(matrix: pandas.DataFrame, bandgap_bins: "list[float]"):
    """Create the bandgap average and bandgap standard deviation for every functional group column (category) in the given matrix"""

    averages: list[float] = [
        float(numpy.average(a=bandgap_bins, weights=matrix[fg])) 
        for fg in matrix.columns
    ]

    variances: list[float] = [
        float(numpy.average(a=[(bandgap-averages[i])**2 for bandgap in bandgap_bins], weights=matrix[fg])) 
        for i, fg in enumerate(matrix.columns)
    ]

    return (averages, [math.sqrt(variance) for variance in variances])

##### Execution Time Variable #####
now = datetime.now()
spinner.next()

##### Data Load #####
bandgap_data: pandas.Series = pandas.read_excel(BANDGAPS_PATH, sheet_name="All")["bandgap"]
fg_data: pandas.DataFrame = pandas.read_excel(MAIN_OUTPUT_PATH, sheet_name="exact_data").drop("AminoAcid", axis=1)
bandgap_fg_data: pandas.DataFrame = pandas.concat([fg_data, bandgap_data], axis=1).set_index("Refcode").drop("SMILES", axis=1)

##### Bandgap Bins #####
max_bandgap: int = max(bandgap_fg_data["bandgap"])
bandgap_bins: "list[float]" = [x*0.5 for x in range(0, 2*math.ceil(max_bandgap))]

##### Bandgap Bin Sorted Molecular Structures #####
bandgap_range_sorted_molecule_dataframes: list[pandas.DataFrame] = [
    bandgap_fg_data.loc[
        (bandgap_fg_data['bandgap'] >= bandgap) & 
        (bandgap_fg_data['bandgap'] < bandgap+0.5)
    ]
    for bandgap in bandgap_bins
]
spinner.next()

##### Structure Count Validity Check #####
molecule_set_sizes: list[int] = [len(molecules) for molecules in bandgap_range_sorted_molecule_dataframes]
assert sum(molecule_set_sizes) == len(bandgap_fg_data)

##### Bandgap Row Organization Lookup #####
bandgap_row_lookup: dict[float, int] = {bg_value: 1+row for row, bg_value in enumerate(bandgap_bins)}

########## Functional Group Independent Counting In Bin Sorted Structure Sets ##########
##### Functional Group Independent Counts Data & Matrix #####
independent_data: list[list[int]] = [
    [len(bandgap_molecules_subtable.loc[bandgap_molecules_subtable[f] > 0]) for f in bandgap_fg_data.columns]
    for bandgap_molecules_subtable in bandgap_range_sorted_molecule_dataframes
]
independent_matrix: pandas.DataFrame = pandas.DataFrame(data=independent_data, columns=bandgap_fg_data.columns, index=bandgap_bins)

##### Independent Matrix Bandgap Means and Standard Deviations Calculations Per Functional Group Column #####
independent_means: list[float]
independent_stds: list[float]
independent_means, independent_stds = createStatisticalMetrics(independent_matrix, bandgap_bins)
independent_matrix.loc["total"] = list(independent_matrix.sum())              # type: ignore
independent_matrix.loc["mean"] = list(independent_means)                      # type: ignore
independent_matrix.loc["std"] = list(independent_stds)                        # type: ignore
spinner.next()

##### Independent Matrix Excel Sheet #####
independent_matrix.to_excel(writer, sheet_name="independent_matrix", freeze_panes=(1, 1))
independent_sheet = writer.sheets["independent_matrix"]
independent_sheet.set_column(0, 0, 7)
independent_sheet_columns: list[str] = [str(col) for col in independent_matrix.columns][1:]
for i, col in enumerate(independent_sheet_columns):
    independent_sheet.set_column(i+1, i+1, len(col)+7)
spinner.next()

########## Functional Group Combinational Counting In Bin Sorted Structure Sets ##########
##### Combinational Sets Loop #####
for i, fg_set in enumerate(SETS):

    ##### Combinational Queries #####
    combinational_queries: list[str] = generateCombinations(fg_set)

    ##### Combinational Category Names #####
    combinational_category_names: list[str] = [' '.join(fgs) if len(fgs) > 1 else fgs[0] for fgs in powerset(fg_set) if fgs]

    ##### Combinational Group Counts Data & Matrix #####
    combinational_data: list[list[int]] = [
        [len(bandgap_molecules_subtable.query(q)) for q in combinational_queries]
        for bandgap_molecules_subtable in bandgap_range_sorted_molecule_dataframes
    ]
    combinational_matrix: pandas.DataFrame = pandas.DataFrame(data=combinational_data, columns=combinational_category_names, index=bandgap_bins)
    combinational_matrix["total"] = [sum(d) for d in combinational_data]
    combinational_matrix = combinational_matrix.loc[:, (combinational_matrix != 0).any(axis=0)]
    spinner.next()

    ##### Combinational Matrix Bandgap Means and Standard Deviations Calculations Per Functional Group Relation #####
    combinational_means: list[float]
    combinational_stds: list[float]
    combinational_means, combinational_stds = createStatisticalMetrics(combinational_matrix, bandgap_bins)
    combinational_matrix.loc["total"] = list(combinational_matrix.sum())              # type: ignore
    combinational_matrix.loc["mean"] = list(combinational_means)                      # type: ignore
    combinational_matrix.loc["std"] = list(combinational_stds)                        # type: ignore
    spinner.next()

    ##### Combinational Matrix Excel Sheet #####
    combinations_sheet_name: str = "combinations_matrix_" + str(i)
    combinational_matrix.to_excel(writer, sheet_name=combinations_sheet_name, freeze_panes=(1, 1))
    combinations_sheet = writer.sheets[combinations_sheet_name]
    combinations_sheet.set_column(0, 0, 7)
    combinations_sheet_columns: list[str] = [str(col) for col in combinational_matrix.columns][1:]
    for j, col in enumerate(combinations_sheet_columns):
        combinations_sheet.set_column(j+1, j+1, len(col)+7)
    spinner.next()

    ########## Chart Creation ##########
    ##### Combinational Matrix Chart Setup #####
    stacked_chart = writer.book.add_chart({'type': 'column', 'subtype': 'stacked'})             # type: ignore
    combinational_rows = combinational_matrix.loc[combinational_matrix["total"] > 0]
    combinational_bandgap_bins = [index for index in combinational_rows.index if index not in ["total", "mean", "std"]]
    fg_groups = [group for group in combinational_matrix.columns if group != "total"]
    spinner.next()

    ##### Bandgap Bin Sorted Combinational Functional Group Percentage Label Creation #####
    custom_labels = {fg: [] for fg in fg_groups}
    for index, item in combinational_rows.iterrows():

        ##### Skip Non-Functional Group Categories #####
        if index in ["total", "mean", "std"]: continue

        ##### Bandgap Categorized Combinational Group Percentages Calculation And Labeling Setup #####
        row: list[float] = list(item)
        total: float = row[-1]
        counts: list[float] = row[:-1]
        percentages: list[str] = [str(round((count/total)*100, 2)) + "%" if count else " "  for count in counts]

        ##### Combinational Groups Labeling To Chart Setup #####
        for j, p in enumerate(percentages):
            custom_labels[fg_groups[j]].append({'value': p, 'font': {'color': COLORS[j].label_color}})

    spinner.next()

    ##### Stacked Chart Plot Insertion #####
    start_row = bandgap_row_lookup[combinational_bandgap_bins[0]]
    for j, (name, item) in enumerate(combinational_rows.items()):
        if name == "total": continue
        stacked_chart.add_series({
            'values':       [combinations_sheet_name, start_row, j+1, start_row+len(combinational_bandgap_bins)-1, j+1], 
            'categories':   [combinations_sheet_name, start_row, 0, start_row+len(combinational_bandgap_bins)-1, 0],
            'data_labels':  {'value': True, 'custom': custom_labels[str(name)]},
            'fill':         {'color': COLORS[j].bar_color},
            'name':         str(name),
        })

    ##### Stacked Chart Axis Labeling #####
    stacked_chart.set_x_axis({
        'name': 'Optical Band Gaps (eV)',
        'name_font': {'size': 30, 'bold': True},
        'num_font':  {'italic': True },
    })

    stacked_chart.set_y_axis({
        'name': 'Structure Count',
        'name_font': {'size': 30, 'bold': True},
        'num_font':  {'italic': True },
    })


    ##### Stacked Chart Size Formatting #####
    stacked_chart.set_size({
        'width': 1860,
        'height': 720
    })

    ##### Stacked Chart Title #####
    stacked_chart.set_title({
        'name': combinations_sheet_name
    })

    ##### Stacked Chart Legend #####
    stacked_chart.set_legend({'font': {'size': 13, 'bold': True}})

    ##### Stacked Chart Sheet Insertion #####
    combinations_sheet.insert_chart("A" + str(len(combinational_matrix.index) + 5), stacked_chart)
    spinner.next()

##### Excel Sheet Save ####
writer.close()
spinner.next()

##### Execution Time Evaluation #####
print("\nexecution time = ", datetime.now() - now, "s")
print(ANALYSIS_OUTPUT_PATH, "file created")

