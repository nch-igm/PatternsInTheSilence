#!/usr/bin/env python

"""
<Purpose>
 
Usage: python my_script.py <in_filename> <out_filename>
Reads: in_filename
Writes: out_filename
 
Usage example:
<command line execution>
 
Modifications:
 <date> - <mod>
"""

import argparse
import pandas as pd
import importlib.util
from generate_matched_context_distribution import get_matched_sampling

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

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
    parser.add_argument('--target-size', '-s', help="Set if group size is known (default is to take minimum group size)",
                        type=int,
                        required=False)
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/1_mutual_information/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    #Read in distribution
    context_df = pd.read_csv(args.codon_file,
                             delimiter="\t",
                             dtype={"CHR":str},
                             usecols=["CHR", "POS", "REF", "REF_Codon", 
                                      "REF_Sequence", "REF_AminoAcid", "REF_AminoAcid_sub"])

    
    #List of target groups
    if args.target_list == ["codons_with_sub_syn_noSTOP"] :
        target_list = ccv.codons_with_sub_syn_noSTOP
    elif args.target_list == ["aa_sub_with_syn_noSTOP"] :
        target_list = ccv.aa_sub_with_syn_noSTOP
    else :
        target_list = args.target_list
    
    #If pre-set target group is listed use that
    if args.target_size is not None :
        target_size = args.target_size
        print("Using set group size:", target_size)
    else :
        #Take original distribution (subset to target list)    
        group_sizes = context_df.value_counts(args.target_column).loc[target_list]
        #-select minimum size
        target_size = group_sizes.min()
        print("Starting group sizes:", group_sizes)
        print("Using minimum group size:", target_size)

    target_group_sizes = pd.Series(target_size,
                                   index=target_list)
    
    print("Size of target data frame:", target_size*len(target_list))
    
    context_sampled_df = get_matched_sampling(context_df,
                                              args.target_column,
                                              target_group_sizes)
    
    print(context_sampled_df.value_counts(args.target_column))
    
    del context_df
    
    sampled_filename = (args.output_folder+
                        args.output_tag+".tsv")
    
    context_sampled_df.to_csv(sampled_filename,
                              sep="\t",
                              index=False)
    
    
    
