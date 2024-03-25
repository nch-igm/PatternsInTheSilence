#!/usr/bin/env python

"""

Sub-sample source codon table to match distribution of codon records 
in the reference codon table
 
Usage: python generate_matched_context_distribution.py \
    -f <reference codon table> -l <aa list> \
        -a <aa column> -f2 <source codon table> \
            -o <output folder> -t <output tag>
Reads: reference codon_tale, source codon table
Writes: <output folder>/<output file>.tsv

"""

import argparse
import pandas as pd
import importlib.util

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

def get_matched_sampling (df, column, group_sizes) :
    
    """
    sample df based on values of column into sizes specified in group_sizes
    """
    
    sampled_df_groups = []
    for group_label, df_group in df.groupby(column) :
        if group_label not in group_sizes :
            continue
        else :
            try:
                sampled_df_group = df_group.sample(n=group_sizes.loc[group_label],
                                                replace=False,
                                                random_state=10)
            except ValueError :
                print(group_label, "want:", group_sizes.loc[group_label], 
                    "- only have:", df_group.value_counts(column).loc[group_label])
                exit(1)
            sampled_df_groups.append(sampled_df_group)
        
    sampled_df = pd.concat(sampled_df_groups, axis=0)
    return sampled_df

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--codon-file', '-f', help="Path to TSV containing Codon and Sequence Context information",
                        required=True)
    parser.add_argument('--target-list', '-l', help="Flag which amino acid list to use (use global variable module label), or enter list",
                        nargs="+",
                        default="codons_with_sub_syn_noSTOP")
    parser.add_argument('--target-column', '-a', help="Name of column in dataframe containing the selected amino acid labels",
                        type=str,
                        default="REF_Codon")
    parser.add_argument('--codon-file2', '-f2', help="Path to TSV containing Codon and Sequence Context information to be sampled from",
                        required=True)
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/1_mutual_information/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    context_df = pd.read_csv(args.codon_file,
                             delimiter="\t",
                             dtype={"CHR":str},
                             usecols=["CHR", "POS", "REF", "REF_Codon", 
                                      "REF_Sequence", "REF_AminoAcid", "REF_AminoAcid_sub"])
    
    if args.target_list == ["codons_with_sub_syn_noSTOP"] :
        target_list = ccv.codons_with_sub_syn_noSTOP
    else :
        target_list = args.target_list
        
    target_group_sizes = context_df.value_counts(args.target_column).loc[target_list]
    
    print("Size of target data frame:", context_df.shape[0])
    print("Target distribution:", target_group_sizes)
    del context_df
    
    context_df2 = pd.read_csv(args.codon_file2,
                             delimiter="\t",
                             dtype={"CHR":str},)
    
    context2_sampled_df = get_matched_sampling(context_df2,
                                               args.target_column,
                                               target_group_sizes)
    
    print(context2_sampled_df.value_counts(args.target_column))
    
    del context_df2
    
    sampled_filename = (args.output_folder+
                        args.output_tag+".tsv")
    
    context2_sampled_df.to_csv(sampled_filename,
                               sep="\t",
                               index=False)
    
    
    
