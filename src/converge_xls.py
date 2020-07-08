import os
#import csv
#from xlsxwriter.workbook import Workbook
import sys
import pandas as pd
#counts_file = sys.argv[1]
#normalized = sys.argv[2]

print("arg list = "+str(sys.argv))

"""
workbook = Workbook("combined_counts.xlsx")
for tsv_file in sys.argv[1:]:
    worksheet = workbook.add_worksheet()
    with open(tsv_file, 'rt') as f:
        reader = csv.reader(f)
        for r, row in enumerate(reader):
            for c, col in enumerate(row):
                worksheet.write(r, c, col)
    workbook.close()
"""


writer = pd.ExcelWriter('readcount_report.xlsx') # Arbitrary output name
for csvfilename in sys.argv[1:]:
    df = pd.read_csv(csvfilename, sep='\t')
    #name = csvfilename.split('/')[-1]
    name = csvfilename
    df.to_excel(writer,sheet_name=name)
writer.save()
