#!/usr/bin/env python

"""

Read in two sets of variant files and join to one file
 
Usage: python combine_scored_variant_files.py \
    -f1 <search string1> -f2 <search string2> \
        -l <aa list> -o <output path>
        
Reads: files matching <search string 1>, files matching <search string 2>
Writes: <output path>
 
"""

import pandas as pd
import glob
import argparse

import importlib.util

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

def combine_on_aa_files (search_string_1, search_string_2, aa_list,
                         merge_keys = ["CHR", "POS", "REF", "ALT"]) :
    
    """
    get filenames for search_string_1, search_string_2, keep 
    records on aa_list, merge each table on merge_keys
    """

    files1 = glob.glob(search_string_1)
    files2 = glob.glob(search_string_2)

    print("Found files in first search:", files1)
    print("Found files in second search:", files2)
    
    merged_dfs = []
    
    for amin in aa_list :
        print("Reading in files for Amino Acid", amin)
        file1 = [x for x in files1
                 if "AminoAcid"+amin in x]
        file2 = [y for y in files2
                 if "AminoAcid"+amin in y]

        if len(file1) > 1:
            print(" Can't uniquely specify file:", files1)
        else :
            pass
        if len(file2) > 1:
            print(" Can't uniquely specify file:", files2)
        else :
            pass

        df1 = pd.read_csv(file1[0],
                          sep="\t",
                          dtype={"CHR":str})
        df2 = pd.read_csv(file2[0],
                          sep="\t",
                          dtype={"CHR":str})

        df_merged = df1.merge(df2,
                              on=merge_keys,
                              how="inner")

        print("Checking key matching:")
        print(" df1:", df1.shape)
        print(" df2:", df2.shape)
        print(" merged:", df_merged.shape)

        merged_dfs.append(df_merged)

    all_df_merged = pd.concat(merged_dfs, ignore_index=True)
    print(all_df_merged.shape)

    return all_df_merged

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--file-source-1', '-f1', help="Search string to grab desired per-Amino Acid files",
                        type=str,
                        required=True)
    parser.add_argument('--file-source-2', '-f2', help="Search string to grab desired per-Amino Acid files",
                        type=str,
                        required=True)
    parser.add_argument('--aa-list', '-l', help="Flag which amino acid list to use (all, at_cp3), or enter list",
                        nargs="+",
                        default="all")
    parser.add_argument('--output-file', '-o', help="Path/filename to label output files with",
                        type=str,
                        default="./")
    args = parser.parse_args()

    if args.aa_list == ["all"] :
        aa_list = ccv.aa_with_syn_noSTOP
    elif args.aa_list == ["at_cp3"] :
        aa_list = ccv.aa_sub_with_syn_noSTOP
    else :
        aa_list = args.aa_list

    #Read in variant file
    merged_df = combine_on_aa_files(args.file_source_1,
                                    args.file_source_2,
                                    aa_list)

    print("Merged columns:", merged_df.columns)

    merged_df.to_csv(args.output_file,
                     sep="\t",
                     index=False)
