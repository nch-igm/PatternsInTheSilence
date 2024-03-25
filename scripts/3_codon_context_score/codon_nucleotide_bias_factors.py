#!/usr/bin/env python

""" 

From codon record table, get counts for codons and context nucleotides,
calculate bias factors for each combination of central codon, context position,
and nucleotide

Usage: python codon_nucleotide_bias_factors.py -f <codon_record_tsv> \
    -l <aa list> -a <aa column in file> -p <# of positions in context> \
        -pr <position range> -o <output folder> \
            -t <output tag>
Read: codon_record_tsv
Write: <output folder>/codon_nuc_mi_bias_factors_<aa column in file-REF>_<position range>nt_<output tag>.tsv

"""

import numpy as np
import pandas as pd
import argparse

import importlib.util
import sys

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

def calc_codon_nuc_counts_tables (codon_df,
                                  amino_acid_list,
                                  amino_acid_col,
                                  num_pos) :
    """
    from codon record table, construct table with counts of each combination of
    central codon, nucleotide and positions
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
    aa_codon_count_nx_pi_dict = {}
    
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
        aa_codon_count_nx_pi_dict[amin] = count_codon_nx_pi
        
    return aa_codon_count_nx_pi_dict

def calc_codon_nuc_mi_factors (aa_codon_count_nx_pi_dict,
                               third_base_pos,
                               pos_range):
    """
    calculate bias factors on given range from count table
    """

    #--Storage
    codon_nuc_bias_rows = []

    #--Indices
    pos_array = (list(range(third_base_pos-3-pos_range,
                            third_base_pos-3))+
                 list(range(third_base_pos,
                            third_base_pos+pos_range)))
    pos_labels = (["N"+str(x) for x in range(-1*pos_range,0)]+
                  ["N+"+str(x) for x in range(1,pos_range+1)])
    
    for a_i, amin in enumerate(aa_codon_count_nx_pi_dict) :
        print("Amino acid:", amin)
        #how many codons = AA
        syn_codons = ccv.aa_comb_syn_codons[amin]
        num_syn_codons = len(syn_codons)

        count_codon_nx_pi = aa_codon_count_nx_pi_dict[amin]
        #Calculate marginal counts
        count_aa = np.sum(count_codon_nx_pi[:,:,0])
        #-N(N_x) per position, per base
        #-shape base x pos
        count_nx_pi = np.sum(count_codon_nx_pi, axis=0)
        #-N(C)
        #-shape codon x 1
        count_codon = np.sum(count_codon_nx_pi[:,:,0], axis=1).reshape(num_syn_codons,1)
        #-N(C,N_x) per position
        #-shape codon x base x pos
        #-count_codon_nx_pi
        
        #Calculate probabilities
        p_codon = count_codon/count_aa
        p_nx_pi = count_nx_pi/count_aa
        p_codon_nx_pi = count_codon_nx_pi/count_aa

        # loop through codons corresponding to this aa
        for c_i, codon in enumerate(syn_codons) :
            # loop through positions in sequence context
            for p_i, pos in enumerate(pos_array) :
                #loop through each base
                for base in ccv.base_dict :
                    codon_nuc_pos_count = count_codon_nx_pi[c_i, ccv.base_dict[base], pos]

                    p_c = p_codon[c_i][0]
                    p_x_n = p_nx_pi[ccv.base_dict[base], pos]
                    p_c_x_n = p_codon_nx_pi[c_i, ccv.base_dict[base], pos]
                    
                    if codon_nuc_pos_count > 0 :
                        mi_factor = p_c_x_n*(np.log(p_c_x_n)
                                               -np.log(p_c)-np.log(p_x_n))
                    else :
                        mi_factor = 0
                        
                    codon_nuc_bias_row = {"AminoAcid": amin,
                                          "AminoAcid_count": int(count_aa),
                                          "Codon": codon,
                                          "Position": pos_labels[p_i],
                                          "Nucleotide": base,
                                          "P_Codon": p_c,
                                          "Pa_Pos_N": p_x_n,
                                          "Pa_Codon_Pos_N": p_c_x_n,
                                          "bias_factor": mi_factor}
                    codon_nuc_bias_rows.append(codon_nuc_bias_row)

    codon_nuc_bias_df = pd.DataFrame(codon_nuc_bias_rows)
    
    return codon_nuc_bias_df

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
    parser.add_argument('--pos-range', '-pr', help="Number of flanking nucleotides to calculate factors on",
                        type=int,
                        default=12)
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/3_codon_context_score/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()

    context_df = pd.read_csv(args.codon_file,
                             delimiter="\t",
                             dtype={"CHR":str})

    if args.aa_list == ["all"] :
        aa_list = ccv.aa_with_syn_noSTOP
    elif args.aa_list == ["at_cp3"] :
        aa_list = ccv.aa_sub_with_syn_noSTOP
    else :
        aa_list = args.aa_list

    print("Processing these amino acids:", aa_list)

    #Calculate count tables
    codon_nuc_counts_dict = calc_codon_nuc_counts_tables (context_df,
                                                          aa_list,
                                                          args.aa_column,
                                                          args.num_pos)

    #Calculate bias factors
    third_base_pos = int((args.num_pos-1)/2+1)
    codon_nuc_bias_df = calc_codon_nuc_mi_factors (codon_nuc_counts_dict,
                                                   third_base_pos,
                                                   args.pos_range)

    #Save dataframe
    print("Saving data frame.")
    factor_filename = (args.output_folder+"codon_nuc_mi_bias_factors_"+
                       args.aa_column[4:]+"_"+
                       str(args.pos_range)+"nt"+
                       args.output_tag+".tsv")
    codon_nuc_bias_df.to_csv(factor_filename,
                             sep="\t", index=False)
