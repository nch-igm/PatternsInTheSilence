#!/usr/bin/env python

"""

Generate shuffles of the codon context table and measure codon-nucleotide mutual information
for each shuffled version of the table
 
Usage: python mutual_information_on_shuffled_contexts.py -f <codon_record_tsv> \
    -l <aa list> -a <aa column in file> -p <# of positions in context> \
        -ns <number of shuffles> --seed-start <seed> \
            -o <output folder> -t <output tag>
Reads: codon_record_tsv
Writes: <output folder>

"""

import argparse
import pandas as pd
import numpy as np
from mutual_information_codon_nuc import calc_mutual_information_tables_codon_nuc
from sklearn.utils import shuffle
import importlib.util

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--codon-file', '-f', help="Path to TSV containing Codon and Sequence Context information",
                        required=True)
    parser.add_argument('--aa-list', '-l', help="Flag which amino acid list to use (all, at_cp3), or enter list",
                        nargs="+",
                        default="all")
    parser.add_argument('--aa-column', '-a', help="Name of column in dataframe containing the selected amino acid labels",
                        type=str,
                        default="REF_AminoAcid")
    parser.add_argument('--num-pos', '-p', help="Number of positions in sequence context field",
                        type=int,
                        default=101)
    parser.add_argument('--num-shuffles', '-ns', help="Number of times to shuffle contexts",
                        type=int,
                        default=1000)
    parser.add_argument('--seed-start',type=int)
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/1_mutual_information/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    context_df = pd.read_csv(args.codon_file,
                             delimiter="\t",
                             dtype={"CHR":str},
                             usecols=["CHR", "POS", "REF", "NM_ID",
                                      "REF_Codon", "REF_Sequence", "REF_AminoAcid", "REF_AminoAcid_sub"])

    if args.aa_list == ["all"] :
        aa_list = ccv.aa_with_syn_noSTOP
    elif args.aa_list == ["at_cp3"] :
        aa_list = ccv.aa_sub_with_syn_noSTOP
    else :
        aa_list = args.aa_list
        
    context_xAA = context_df.groupby(args.aa_column)
    random_seeds = np.arange(args.seed_start, (args.seed_start+args.num_shuffles))
    
    #loop through shuffles
    for i in range(args.num_shuffles) :
        print("Running shuffle", i, "with seed", random_seeds[i])

        context_shuffle_df = context_df.copy()
        #loop through amino acid list
        for amin in aa_list :
            aa_idx = context_xAA.groups[amin]
            
            #shuffle contexts among records for that amino acid
            context_shuffle_df.loc[aa_idx, "REF_Sequence"] = \
                shuffle(context_shuffle_df.loc[aa_idx]["REF_Sequence"].values,
                        random_state=random_seeds[i])
            
        #with shuffled data set, calculate mutual information
        a_df, b_df, c_df, mi_codon_nuc_pos_shuffle_df = \
            calc_mutual_information_tables_codon_nuc(context_shuffle_df,
                                                     aa_list,
                                                     args.aa_column,
                                                     args.num_pos)
        #save full set of mutual information results for 
        # this shuffle of the data set
        mi_shuffle_filename = (args.output_folder+"mut_info_codon_nuc_pos_"+
                               str(args.num_pos)+"bp"+
                               args.output_tag+"_"+
                               "shuffle"+str(i)+".tsv")
        mi_codon_nuc_pos_shuffle_df.to_csv(mi_shuffle_filename,
                                           sep="\t", index=True)


