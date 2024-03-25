#!/usr/env python

"""

Take GFF3 file with exon entries and summarize the genic bounds
 
Usage: python annotate_transcripts_with_exons_genomic_range.py --gff-file <GFF file> --file-root <out_root>
Reads: GFF file
Writes: <out_root>_Trans_genomic_coords.tsv, <out_root>_Exon_pos.tsv

"""

import argparse
import pandas as pd
import numpy as np
import sys

pd.set_option("display.max_columns", None)

def get_exon_info (x) :
    """
    parse GFF info fields
    """
    fields = x.split(";")
    field_dict = {y.split("=")[0]:y.split("=")[1] for y in fields}

    Target_fields = field_dict["Target"].split(" ")
    NM_ID = Target_fields[0]
    Trans_Start = int(Target_fields[1])
    Trans_Stop = int(Target_fields[2])
    
    num_mismatch = int(field_dict["num_mismatch"])
    gap_count = int(field_dict["gap_count"])
    num_ident = int(field_dict["num_ident"])

    info_dict = {"NM_ID": NM_ID,
                 "Exon_Trans_Start": Trans_Start,
                 "Exon_Trans_Stop": Trans_Stop,
                 "num_mismatch": num_mismatch,
                 "gap_count": gap_count,
                 "num_ident": num_ident}
    
    return info_dict

def get_key (r, key_tuple = ["NM_ID", "UCSC_CHR", "Sense"]) :
    """
    return a key to identify transcripts
    """
    row_key = tuple([r[k] for k in key_tuple])

    return row_key

def get_transcript_summary (gcf_df) :
    """
    go through GFF file and summarize transcripts from exon records that map to 
    the transcript ID
    """
    #In the GFF file, the exon entries are listed in groups, and sorted
    # by transcript position
    #Exon_Start < Exon_Stop AND Exon_Trans_Start < Exon_Trans_Stop
    # but exons are listed in transcript-order
    
    transcript_dict = {}
    transcript_exons_dict = {}
    last_exon_key = (None) #store which transcript/loc was processed last
    last_exon_trans_stop = 1000
    exon_counter = 1

    for idx, row in gcf_df.iterrows() : #iterate through the exon rows
        this_exon_key = get_key(row) #NM_ID, UCSC_CHR, Sense
        #starts and stops are always listed in increasing order
        this_exon_trans_start = row["Exon_Trans_Start"]
        
        if (this_exon_key != last_exon_key and
            this_exon_trans_start <= last_exon_trans_stop):
            #start new transcript
            #-get genomic coordinate for 5' end of mapped transcript
            if row["Sense"] == "+" :
                trans_pos_start = int(row["Exon_Start"])
            else :
                trans_pos_start = int(row["Exon_Stop"])

            this_transcript_dict = {"NM_ID":row["NM_ID"],
                                    "UCSC_CHR":row["UCSC_CHR"],
                                    "Sense":row["Sense"],
                                    "POS_Start":trans_pos_start,
                                    "Exon_Count":1,
                                    "POS_Left":row["Exon_Start"],
                                    "POS_Right":row["Exon_Stop"],
                                    "Total_Exonic_Length_Mapped":row["Exon_Length_Mapped"],
                                    "num_ident":row["num_ident"]}

            #Add to transcript_dict
            this_transcript_key = tuple(list(this_exon_key)+
                                        [trans_pos_start])
            print(">the new transcript key is", this_transcript_key)
            
            transcript_dict[this_transcript_key] = this_transcript_dict
            last_exon_key = this_exon_key ##
            last_exon_trans_stop = row["Exon_Trans_Stop"]
            
            gcf_df.loc[idx, "Exon_Number"] = 1
            exon_counter = 2 
        else :
            #this is the next exon in the current transcript
            transcript_dict[this_transcript_key]["Exon_Count"] += 1
            if row["Sense"] == "+" :
                transcript_dict[this_transcript_key]["POS_Right"] = row["Exon_Stop"]
            else :
                transcript_dict[this_transcript_key]["POS_Left"] = row["Exon_Start"]
            transcript_dict[this_transcript_key]["Total_Exonic_Length_Mapped"] += row["Exon_Length_Mapped"]
            last_exon_trans_stop = row["Exon_Trans_Stop"]

            gcf_df.loc[idx, "Exon_Number"] = exon_counter
            exon_counter += 1
            
    gcf_transcript_df = pd.DataFrame.from_dict(transcript_dict,
                                               orient="index")
    print("#Transcript summary:")
    print(gcf_transcript_df.head())
    print("#Updated exon summary:")
    print(gcf_df.head())

    return gcf_transcript_df, gcf_df

if __name__ == "__main__" :
    #Get inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff-file", "-g", help="Path to GFF3 file with exon entries")
    parser.add_argument("--file-root", "-f", help="Path/label to apply to start of output files")
    args = parser.parse_args()
    
    #Read in file
    gcf_df = pd.read_csv(args.gff_file,
                         sep="\t",
                         comment="#",
                         header=None,
                         names=["UCSC_CHR", "Database", "AlignmentType",
                                "Exon_Start", "Exon_Stop", "Exon_Length_Mapped",
                                "Sense", "Qual", "Info"])

    #Extract and add Target information
    gcf_df["Info_dict"] = gcf_df["Info"].apply(get_exon_info)
    gcf_exon_df = pd.concat([gcf_df,
                             gcf_df["Info_dict"].apply(pd.Series)],
                            axis=1)

    print(gcf_exon_df.shape)
    print(gcf_exon_df.head())

    #Summarize some of the transcript-level attributes
    gcf_exon_df["Exon_Number"] = 0
    gcf_transcript_df, gcf_exon_df = get_transcript_summary(gcf_exon_df)

    #Save files
    gcf_transcript_output_filename = args.file_root+"_Trans_genomic_coords.tsv"
    gcf_exon_output_filename = args.file_root+"_Exon_pos.tsv"

    gcf_exon_df.\
        drop(columns=["Database", "AlignmentType", "Qual", "Info_dict"]).\
        to_csv(gcf_exon_output_filename,
               sep="\t",
               index=False)

    gcf_transcript_df.to_csv(gcf_transcript_output_filename,
                             sep="\t",
                             index=False)
    
