#!/usr/bin/env python

""" 

Calculate mutual information between the central codon and the surrounding nucleotides 
and output per-Amino Acid tables for mutual information and the composite shannon entropies

Usage: python mutual_information_codon_nuc.py -f <codon_record_tsv> -l <aa list> \
    -a <aa column in file> -p <# of positions in context> \
        -o <output folder> -t <output tag>
Reads: codon_record_tsv
Writes: <output folder>shannon_entropy_codon<output tag>.tsv
    <output folder>shannon_entropy_nuc_pos_<num pos>bp<output tag>.tsv
    <output folder>shannon_entropy_codon_nuc_pos_<num pos>bp<output tag>.tsv
    <output folder>mut_info_codon_nuc_pos_<num pos>bp<output tag>.tsv
    
"""

import numpy as np
import pandas as pd
import argparse

import importlib.util
import sys

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

def calc_mutual_information_tables_codon_nuc (codon_df,
                                              amino_acid_list,
                                              amino_acid_col,
                                              num_pos) :
    """
    calculate mutual information between codon and nucleotide distributions
    using REF_Sequence column in codon_df for specified set of amino acids
    """
    #---Which AA and codons ---
    num_aa = len(amino_acid_list)
    codons_in_list = [x for x in ccv.codon_aa_dict if 
                      ((ccv.codon_aa_dict[x] in amino_acid_list) or
                       (ccv.codon_aa_sub_dict[x] in amino_acid_list))]
    codons_in_list_indx_dict = {x:i for i,x in enumerate(codons_in_list)}
    num_codons = len(codons_in_list)
    num_bases = ccv.num_bases
    
    #---Storage arrays---
    sh_entropy_codon_list = []
    sh_entropy_nuc_pos_list = []
    sh_entropy_codon_nuc_pos_list = []
    
    #---Group records Amino Acid x Codon---
    codon_xAAxCodon = codon_df.groupby([amino_acid_col, "REF_Codon"])

    #---Calculate counts----
    # split reference sequence in array of letters
    codon_df["REF_Sequence_array"] = codon_df["REF_Sequence"].apply(lambda x: list(x))

    # loop through amino acids in list
    for a_i, amin in enumerate(amino_acid_list) :
        print("Amino acid:", amin)
        #how many codons = AA
        syn_codons = ccv.aa_comb_syn_codons[amin]
        num_syn_codons = len(syn_codons)
        #set up storage
        #- count:
        # Codon x Base observed x At position
        count_codon_nx_pi = np.zeros((num_syn_codons, 
                                      num_bases,
                                      num_pos))
        # loop through codons corresponding to this aa
        for c_i, codon in enumerate(syn_codons) :
            #get the index for this AA, Codon combination
            loc_aa_codon = (amin, codon)
            #get the indexes of all matching codon records
            group_aa_codon = codon_xAAxCodon.groups[loc_aa_codon]
            #stack the context arrays
            array_codon_pos = np.vstack(codon_df.loc[group_aa_codon, "REF_Sequence_array"])
            #convert array to data frame
            df_codon_pos = pd.DataFrame(array_codon_pos)
            # loop through positions in sequence context
            for p_i, pos in enumerate(df_codon_pos.columns) :
                #for each (position) column, count the number of observations
                # for each base
                base_pos_counts = df_codon_pos[pos].value_counts()
                #loop through each base
                for base in base_pos_counts.index :
                    #store in overall array
                    count_codon_nx_pi[c_i, ccv.base_dict[base], p_i] = base_pos_counts[base]
                    #Calculate probabilities
                    #-N(AA)
                    #-shape aa x 1
        count_aa = np.sum(count_codon_nx_pi[:,:,0])
        #-N(N_x) per position, per base
        #-shape pos x base
        count_nx_pi = np.sum(count_codon_nx_pi, axis=0)
        #-N(C)
        #-shape codon x 1
        count_codon = np.sum(count_codon_nx_pi[:,:,0], axis=1).reshape(num_syn_codons,1)
        #-N(C,N_x) per position
        #-shape codon x base x pos
        #-count_codon_nx_pi
        
        #Calculate shannon entropies
        p_codon = count_codon/count_aa
        #(1,)
        h_codon = -1.*np.sum(np.multiply(p_codon,
                                         np.log(p_codon)))
        p_nx_pi = count_nx_pi/count_aa
        #(window_size,)
        h_nx_pi = -1.*np.sum(np.multiply(p_nx_pi,
                                         np.log(p_nx_pi)), axis=0)
        p_codon_nx_pi = count_codon_nx_pi/count_aa
        #(window_size,)
        h_codon_nx_pi = -1.*np.sum(np.multiply(p_codon_nx_pi,
                                               np.log(p_codon_nx_pi)), axis=(0,1))
    
        sh_entropy_codon_list.append(h_codon)
        sh_entropy_nuc_pos_list.append(h_nx_pi)
        sh_entropy_codon_nuc_pos_list.append(h_codon_nx_pi)
        
    sh_entropy_codon_df = pd.DataFrame(np.array(sh_entropy_codon_list), 
                                       columns=["h"],
                                       index=amino_acid_list)
    sh_entropy_nuc_pos_df = pd.DataFrame(np.array(sh_entropy_nuc_pos_list),
                                         columns=["h_p"+str(n) for n in range(num_pos)],
                                         index=amino_acid_list)
    sh_entropy_codon_nuc_pos_df = pd.DataFrame(np.array(sh_entropy_codon_nuc_pos_list),
                                               columns=["h_p"+str(n) for n in range(num_pos)],
                                               index=amino_acid_list)

    #Calculate and store Mutual Information as dataframe
    mi_codon_nuc_pos_array = ((sh_entropy_codon_df.values +
                               sh_entropy_nuc_pos_df.values) -
                              sh_entropy_codon_nuc_pos_df.values)
    center_pos = int((num_pos-1)/2+1)
    mi_codon_nuc_pos_array[:,(center_pos-3):(center_pos)] = 0.0
    mi_codon_nuc_pos_df = pd.DataFrame(mi_codon_nuc_pos_array,
                                       index=sh_entropy_codon_nuc_pos_df.index,
                                       columns=["mi_p"+str(n) for n in range(num_pos)])

    
    return sh_entropy_codon_df, sh_entropy_nuc_pos_df, sh_entropy_codon_nuc_pos_df, mi_codon_nuc_pos_df

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
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/1_mutual_information/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()

    context_df = pd.read_csv(args.codon_file,
                             delimiter="\t")

    if args.aa_list == ["all"] :
        aa_list = ccv.aa_with_syn_noSTOP
    elif args.aa_list == ["at_cp3"] :
        aa_list = ccv.aa_sub_with_syn_noSTOP
    else :
        aa_list = args.aa_list

    print("Processing these amino acids:", aa_list)

    #Calculate entropy tables
    sh_entropy_codon_df, sh_entropy_nuc_pos_df, sh_entropy_codon_nuc_pos_df, mi_codon_nuc_pos_df = \
        calc_mutual_information_tables_codon_nuc(context_df,
                                                 aa_list,
                                                 args.aa_column,
                                                 args.num_pos)
    
    #Save dataframes
    print("Saving data frames.")
    codon_filename = args.output_folder+"shannon_entropy_codon"+args.output_tag+".tsv"
    nuc_pos_filename = args.output_folder+"shannon_entropy_nuc_pos_"+str(args.num_pos)+"bp"+args.output_tag+".tsv"
    codon_nuc_pos_filename = args.output_folder+"shannon_entropy_codon_nuc_pos_"+str(args.num_pos)+"bp"+args.output_tag+".tsv"
    mi_codon_nuc_pos_filename = args.output_folder+"mut_info_codon_nuc_pos_"+str(args.num_pos)+"bp"+args.output_tag+".tsv"
    
    sh_entropy_codon_df.to_csv(codon_filename,
                               sep="\t", index=True)
    sh_entropy_nuc_pos_df.to_csv(nuc_pos_filename,
                                 sep="\t", index=True)
    sh_entropy_codon_nuc_pos_df.to_csv(codon_nuc_pos_filename,
                                       sep="\t", index=True)
    mi_codon_nuc_pos_df.to_csv(mi_codon_nuc_pos_filename,
                               sep="\t", index=True)
