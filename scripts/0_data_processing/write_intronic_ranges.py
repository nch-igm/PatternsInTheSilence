#!/usr/env python3

"""
Read in table of exon coordinates, write out table of intron coordinates.

usage: python write_intronic_ranges.py <exon_coord_file> <intron_coord_file>
reads: exon_coord_file
writes: intron_coord_file

"""

import sys
import pandas as pd

#Read in file where each row is an exon and exons are printed in
# ascending order on the mRNA
exon_coord_file = sys.argv[1]
intron_coord_file = sys.argv[2]

exon_coord_df = pd.read_csv(exon_coord_file,
                            sep="\t")

current_intron_num = 0
intron_rows = []
last_row = 0
for idx, row in exon_coord_df.iterrows() : #iterate through the exon rows
    #If we reach new transcript start over
    if (current_intron_num == 0 or
        row["Exon_Number"] == 1) :
        last_row = row
        current_intron_dict = {}
        current_intron_num = 1
    else :
        #Compare to previous intron
        if row["Sense"] == "+" :
            intron_start = last_row["Exon_Stop"]+1
            intron_stop = row["Exon_Start"]-1
        else :
            intron_start = row["Exon_Stop"]+1
            intron_stop = last_row["Exon_Start"]-1

        intron_len = intron_stop-intron_start+1
        current_intron_dict = {"UCSC_CHR": row["UCSC_CHR"],
                               "Intron_Start":intron_start,
                               "Intron_Stop":intron_stop,
                               "Intron_Length_Genomic":intron_len,
                               "Sense":row["Sense"],
                               "NM_ID":row["NM_ID"],
                               "Intron_Number":current_intron_num}
        intron_rows.append(current_intron_dict)

        current_intron_num += 1
        last_row = row

intron_coord_df = pd.DataFrame(intron_rows)
intron_coord_df.to_csv(intron_coord_file,
                       sep="\t", index=False)
