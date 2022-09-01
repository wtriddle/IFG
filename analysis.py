"""Analyzes a functional group dataset and a bandgap dataset for statistically relevant trends"""

import pandas
import math
import numpy
from config import STATS_OUTPUT_PATH
from itertools import chain, combinations
from datetime import datetime

##### User Defined Relational Functional Group Sets #####
SETS = [
    ["PrimaryAmine", "SecondaryAmine", "TertiaryAmine"],
    ["Aromatic Rings", "Non Aromatic Rings"],
]

##### Relational Functional Groups Color Schemes #####
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
writer = pandas.ExcelWriter(STATS_OUTPUT_PATH, engine="xlsxwriter")

##### Powerset Function #####
def powerset(iterable):
    """Produce all subsets, including the full set, of a given iterable"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


##### Relational Query Generator Function #####
def generateQueries(fg_set: "list[str]") -> "list[str]":
    """Generate all relational functional group queries for a pandas dataframe using the functional group names listed in fg_set"""
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
def createStatistics(matrix: pandas.DataFrame, bandgap_bins: "list[float]"):
    """Create the bandgap average and bandgap standard deviation for every functional group column (category) in the given matrix"""
    averages: list[float] = [numpy.average(a=bandgap_bins, weights=matrix[fg]) for fg in matrix.columns]
    variances = [numpy.average(a=([(bandgap-averages[i])**2 for bandgap in bandgap_bins]), weights=matrix[fg]) for i, fg in enumerate(matrix.columns)]
    return (averages, [math.sqrt(variance) for variance in variances])


def main():

    ##### Execution Time Variable #####
    now = datetime.now()

    ##### Data Load #####
    bandgap_data: pandas.Series = pandas.read_excel("CrystalData.xlsx", sheet_name="All")["bandgap"]
    fg_data: pandas.DataFrame = pandas.read_excel("output/output.xlsx", sheet_name="exact_data").drop("AminoAcid", axis=1)
    bandgap_fg_data: pandas.DataFrame = pandas.concat([fg_data, bandgap_data], axis=1).set_index("Refcode")

    ##### Bandgap Bins #####
    max_bandgap: int = max(bandgap_fg_data["bandgap"])
    bandgap_bins: "list[float]" = [x*0.5 for x in range(0, 2*math.ceil(max_bandgap))]

    ##### Bandgap Bin Sorted Molecular Structures #####
    bandgap_bin_sorted_structure_sets: list[pandas.DataFrame] = [
        bandgap_fg_data.loc[
            (bandgap_fg_data['bandgap'] >= bandgap) & 
            (bandgap_fg_data['bandgap'] < bandgap+0.5)
        ]
        for bandgap in bandgap_bins
    ]

    ##### Structure Count Validity Check #####
    structure_set_sizes: list[int] = [len(structures) for structures in bandgap_bin_sorted_structure_sets]
    assert sum(structure_set_sizes) == len(bandgap_fg_data)
    
    ##### Bandgap Row Organization Lookup #####
    bandgap_row_lookup: dict[float, int] = {bg_value: 1+row for row, bg_value in enumerate(bandgap_bins)}

    ########## Functional Group Instance Counting In Bin Sorted Structure Sets ##########
    ##### Functional Group Instance Counts Data & Matrix #####
    instance_data: list[list[int]] = [
        [len(structures.loc[structures[f] > 0]) for f in bandgap_fg_data.columns]
        for structures in bandgap_bin_sorted_structure_sets
    ]
    instance_matrix: pandas.DataFrame = pandas.DataFrame(data=instance_data, columns=bandgap_fg_data.columns, index=bandgap_bins)

    ##### Instance Matrix Bandgap Means and Standard Deviations Calculations Per Functional Group Column #####
    instance_means: list[float]
    instance_stds: list[float]
    instance_means, instance_stds = createStatistics(instance_matrix, bandgap_bins)
    instance_matrix.loc["total"] = list(instance_matrix.sum())              # type: ignore
    instance_matrix.loc["mean"] = list(instance_means)                      # type: ignore
    instance_matrix.loc["std"] = list(instance_stds)                        # type: ignore

    ##### Instance Matrix Excel Sheet #####
    instance_matrix.to_excel(writer, sheet_name="all_matrix", freeze_panes=(1, 1))
    instance_sheet = writer.sheets["all_matrix"]
    instance_sheet.set_column(0, 0, 7)
    for i in range(0, len(instance_matrix.columns)):
        instance_sheet.set_column(i+1, i+1, len(instance_matrix.columns[i])+7)

    ########## Functional Group Relational Counting In Bin Sorted Structure Sets ##########
    ##### Relational Sets Loop #####
    for i, fg_set in enumerate(SETS):

        ##### Relational Queries #####
        relational_queries: list[str] = generateQueries(fg_set)

        ##### Relational Category Names #####
        relational_category_names: list[str] = [' '.join(fgs) if len(fgs) > 1 else fgs[0] for fgs in powerset(fg_set) if fgs]

        ##### Relational Group Counts Data & Matrix #####
        relational_data: list[list[int]] = [
            [len(structures.query(q)) for q in relational_queries]
            for structures in bandgap_bin_sorted_structure_sets
        ]
        relational_matrix: pandas.DataFrame = pandas.DataFrame(data=relational_data, columns=relational_category_names, index=bandgap_bins)
        relational_matrix["total"] = [sum(d) for d in relational_data]
        relational_matrix = relational_matrix.loc[:, (relational_matrix != 0).any(axis=0)]

        ##### Relational Matrix Bandgap Means and Standard Deviations Calculations Per Functional Group Relation #####
        relational_means: list[float]
        relational_stds: list[float]
        relational_means, relational_stds = createStatistics(relational_matrix, bandgap_bins)
        relational_matrix.loc["total"] = list(relational_matrix.sum())              # type: ignore
        relational_matrix.loc["mean"] = list(relational_means)                      # type: ignore
        relational_matrix.loc["std"] = list(relational_stds)                        # type: ignore

        ##### Relational Matrix Excel Sheet #####
        relations_sheet_name: str = "relations_matrix_" + str(i)
        relational_matrix.to_excel(writer, sheet_name=relations_sheet_name, freeze_panes=(1, 1))
        relations_sheet = writer.sheets[relations_sheet_name]
        relations_sheet.set_column(0, 0, 7)
        for j in range(0, len(relational_matrix.columns)):
            relations_sheet.set_column(j+1, j+1, len(relational_matrix.columns[j])+7)

        ########## Chart Creation ##########
        ##### Relational Matrix Chart Setup #####
        stacked_chart = writer.book.add_chart({'type': 'column', 'subtype': 'stacked'})
        relational_rows = relational_matrix.loc[relational_matrix["total"] > 0]
        relational_bandgap_bins = [index for index in relational_rows.index if index not in ["total", "mean", "std"]]
        fg_groups = [group for group in relational_matrix.columns if group != "total"]

        ##### Bandgap Bin Sorted Relational Functional Group Percentage Label Creation #####
        custom_labels = {fg: [] for fg in fg_groups}
        for index, item in relational_rows.iterrows():

            ##### Skip Non-Functional Group Categories #####
            if index in ["total", "mean", "std"]: continue

            ##### Bandgap Categorized Relational Group Percentages Calculation And Labeling Setup #####
            row: list[float] = list(item)
            total: float = row[-1]
            counts: list[float] = row[:-1]
            percentages: list[str] = [str(round((count/total)*100, 2)) + "%" if count else " "  for count in counts]

            ##### Relational Groups Labeling To Chart Setup #####
            for j, p in enumerate(percentages):
                custom_labels[fg_groups[j]].append({'value': p, 'font': {'color': COLORS[j].label_color}})

        ##### Stacked Chart Plot Insertion #####
        start_row = bandgap_row_lookup[relational_bandgap_bins[0]]
        for j, (name, item) in enumerate(relational_rows.iteritems()):
            if name == "total": continue
            stacked_chart.add_series({
                'values':       [relations_sheet_name, start_row, j+1, start_row+len(relational_bandgap_bins)-1, j+1], 
                'categories':   [relations_sheet_name, start_row, 0, start_row+len(relational_bandgap_bins)-1, 0],
                'data_labels':  {'value': True, 'custom': custom_labels[str(name)]},
                'fill':         {'color': COLORS[j].bar_color},
                'name':         str(name),
            })

        ##### Stacked Chart Axis Labeling #####
        stacked_chart.set_x_axis({
            'name': 'Bandgap Ranges',
            'name_font': {'size': 18, 'bold': True},
            'num_font':  {'italic': True },
        })

        stacked_chart.set_y_axis({
            'name': 'Structure Count',
            'name_font': {'size': 18, 'bold': True},
            'num_font':  {'italic': True },
        })


        ##### Stacked Chart Size Formatting #####
        stacked_chart.set_size({
            'width': 1860,
            'height': 720
        })

        ##### Stacked Chart Title #####
        stacked_chart.set_title({
            'name': relations_sheet_name
        })

        ##### Stacked Chart Legend #####
        stacked_chart.set_legend({'font': {'size': 13, 'bold': True}})

        ##### Stacked Chart Sheet Insertion #####
        relations_sheet.insert_chart("A" + str(len(relational_matrix.index) + 5), stacked_chart)

    ##### Excel Sheet Save ####
    writer.save()

    ##### Execution Time Evaluation #####
    print("execution time = ", datetime.now() - now, "s")
    print(STATS_OUTPUT_PATH, "file created")



if __name__ == "__main__":
    main()