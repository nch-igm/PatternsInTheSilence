#!/usr/bin/env python

""" 

Read in data frame with RefSeq transcript IDs and add Entrez Gene IDs, if available

Usage: python annotate_transcript_list_with_entrez.py <transcript_filename> <transcript_annotated_filename>
Reads: transcript_filename
Writes: transcript_annotated_filename

"""

import pandas as pd
import mygene
import sys

#Filenames
#input:
transcript_filename = sys.argv[1]
#output:
transcript_annotated_filename = sys.argv[2]

#Read in file with RefSeq IDs
print("Reading:", transcript_filename)
transcript_df = pd.read_csv(transcript_filename, sep="\t")

#Initialize mygene class
mg = mygene.MyGeneInfo()
#-query with RefSeq transcript IDs
transcript_query_df = mg.querymany(transcript_df["NM_ID"].values, species='human', scopes='refseq', as_dataframe=True)

#-check how many labels are left unannotated
print("dimension of input:", transcript_df.shape)
print(" dimension of output:", transcript_query_df.shape)
print(" dimension of structure with no Entrez Gene:", transcript_query_df[transcript_query_df["entrezgene"].isna()].shape)

#Merge query data frame with original data frame
transcript_anno_df = transcript_df.set_index("NM_ID").\
    merge(transcript_query_df,
          left_index=True,
          right_index=True,
          how="left").\
    reset_index().\
    rename(columns={"index":"NM_ID"})
print("transcript_anno_df:")
print(transcript_anno_df.head())

#Choose final set of columns to save
write_columns = ['NM_ID', 'CHR', 'Sense', 'CDS_Start', 'CDS_Stop', 'CDS_Length',
                 'Trans_Length','entrezgene', 'name', 'symbol']

transcript_anno_write_df = transcript_anno_df[write_columns]

print("writing to:", transcript_annotated_filename)
print(transcript_anno_write_df.head())

transcript_anno_write_df.to_csv(transcript_annotated_filename,
                                sep="\t",
                                index=False)
