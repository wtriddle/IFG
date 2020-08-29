import pandas as pd
from main import identifyFunctionalGroups
import os

df1 = pd.read_excel(
    '../output/FunctionalGroups.xlsx')
df2 = pd.read_excel('../output/FunctionalGroupsTest.xlsx')

print(df1)
print(df2)
diff = df1[df1 != df2].dropna(axis=0, how='all').dropna(axis=1, how='all')
print(diff)


# Retrieve functional groups dataframe
# test = identifyFunctionalGroups(True, True, True)

# # Create pandas excel writer with xlsxwriter as engine
# writer = pd.ExcelWriter('FgTesting.xlsx', engine="xlsxwriter")

# # Covert dataframe to XlsxWriter Excel object
# test['allDf'].to_excel(
#     writer,
#     sheet_name="All Functional Groups",
#     na_rep=0,
#     freeze_panes=(1, 1)
# )
# test['preciseDf'].to_excel(
#     writer,
#     sheet_name="Precise Functional Groups",
#     na_rep=0,
#     freeze_panes=(1, 1)
# )

# # XlsxWriter book/sheet objects
# workbook = writer.book
# allSheet = writer.sheets['All Functional Groups']
# preciseSheet = writer.sheets['Precise Functional Groups']

# # Set sheet column widths
# allSheet.set_column(1, len(test['allDf'].columns), 21, None)
# preciseSheet.set_column(1, len(test['preciseDf'].columns), 21, None)

# writer.save()
