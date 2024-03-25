#!/usr/bin/env python

""" 

For variant records in variant_file, use REF_Sequence and REF_Codon
to calculate REF codon context score. Use REF_Sequence and ALT_Codon
to calculate ALT codon context score. Calculate Delta codon context
score

Usage: python annotate_variant_table_codon_score.py \
    -f <variant_file> -ff <bias factor file> \
        -p <# positions in context> -pr <position range> \
            -o <output folder> -t <output tag>
Reads: variant_file, bias factor file
Writes: <output folder>ssnv_codon_context_score_mi_annotated_<position range>nt_<output tag>.tsv

"""

import argparse
import pandas as pd

def get_pos_labels (pos_range) :
    """
    Column labels for positions in the sequence context
    """

    pos_labels = (["N"+str(x) for x in range(-1*pos_range,0)]+
                  ["N+"+str(x) for x in range(1,pos_range+1)])

    return pos_labels

def annotate_variant_context (variant_df,
                              codon_nuc_bias_df,
                              third_base_pos,
                              pos_range) :
    """
    Add position columns, join codon_nuc_bias_df on 
    codon/nucleotide/position entries to annotate with
    bias factors for each sequence context
    """

    pos_array = (list(range(third_base_pos-3-pos_range,
                            third_base_pos-3))+
                 list(range(third_base_pos,
                            third_base_pos+pos_range)))
    pos_labels = get_pos_labels(pos_range)

    #drop extraneous columns
    aa_columns = [y for y in codon_nuc_bias_df.columns
                  if "AminoAcid" in y]
    codon_nuc_bias_df.drop(columns= (aa_columns + 
                                     ["P_Codon", "Pa_Pos_N", "Pa_Codon_Pos_N"]),
                           inplace=True)

    #Generate nucleotide columns and concat
    pos_columns = [variant_df["REF_Sequence"].apply(lambda x:x[pos])
                   for pos in pos_array]
    pos_columns_df = pd.concat(pos_columns,
                               axis=1)
    pos_columns_df.columns = pos_labels
    print("Generating position columns, indexed:", pos_columns_df.columns)
    
    variant_scored_df = pd.concat([variant_df,
                                   pos_columns_df],
                                  axis=1)
    print("Variant table columns, updated:", variant_scored_df.columns)
    
    for p_i, pos in enumerate(pos_array) :
        pos_label = pos_labels[p_i]
        print("Annotating position:", pos_label, "at index", pos)
        
        #Get subset of codon_nuc_bias_df for this position
        #-select for position
        #-rename column
        codon_nuc_bias_REF_df = codon_nuc_bias_df.query("Position == @pos_label").\
            rename(columns={"Codon":"REF_Codon",
                            "bias_factor":pos_label+"_ref_score",
                            "Nucleotide":pos_label}).drop(columns="Position")
        
        codon_nuc_bias_ALT_df = codon_nuc_bias_df.query("Position == @pos_label").\
            rename(columns={"Codon":"ALT_Codon",
                            "bias_factor":pos_label+"_alt_score",
                            "Nucleotide":pos_label}).drop(columns="Position")

        #Join for REF to add pos_label+"_ref_score" column
        variant_scored_df = variant_scored_df.merge(codon_nuc_bias_REF_df,
                                                    on=["REF_Codon", pos_label],
                                                    how="inner")
        #Join for ALT to add pos_label+"_alt_score" column
        variant_scored_df = variant_scored_df.merge(codon_nuc_bias_ALT_df,
                                                    on=["ALT_Codon", pos_label],
                                                    how="inner")

    return variant_scored_df

def calculate_codon_context_score (variant_scored_df,
                                   pos_range) :
    """
    from variant file annotated with sequence waits, sum to get REF and ALT
    codon context scores and take difference to score sSNV
    """

    pos_labels = get_pos_labels(pos_range)
    
    ref_pos_labels = [x+"_ref_score" for x in pos_labels]
    alt_pos_labels = [x+"_alt_score" for x in pos_labels]
    diff_pos_labels = [x+"_diff_score" for x in pos_labels]

    for p_i, pos in enumerate(pos_labels) :
        variant_scored_df[diff_pos_labels[p_i]] = (variant_scored_df[alt_pos_labels[p_i]] -
                                                   variant_scored_df[ref_pos_labels[p_i]])

    variant_scored_df["sum_ref_context_score"] = variant_scored_df[ref_pos_labels].sum(axis=1)
    variant_scored_df["sum_alt_context_score"] = variant_scored_df[alt_pos_labels].sum(axis=1)
    variant_scored_df["diff_sum_context_score"] = variant_scored_df["sum_alt_context_score"] - variant_scored_df["sum_ref_context_score"]

    return variant_scored_df

def clean_up_score_columns (variant_scored_df,
                            pos_range) :
    """
    Remove per-position columns related to score and nucleotide
    """

    #Remove N columns
    pos_labels = get_pos_labels(pos_range)
    #Individual score columns
    ref_pos_labels = [x+"_ref_score" for x in pos_labels]
    alt_pos_labels = [x+"_alt_score" for x in pos_labels]
    diff_pos_labels = [x+"_diff_score" for x in pos_labels]

    remove_columns = (pos_labels+
                      ref_pos_labels+alt_pos_labels+
                      diff_pos_labels)
    print("Removing following columns from 'clean' version:", remove_columns)
    variant_scored_clean_df = variant_scored_df.drop(columns=remove_columns)

    return variant_scored_clean_df

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--variant-file', '-f', help="Path to TSV containing sSNV and Sequence Context information",
                        required=True)
    parser.add_argument('--factor-file', '-ff', help="Path to TSV containing nucleotide bias factor file",
                        required=True)
    parser.add_argument('--num-pos', '-p', help="Number of positions in sequence context field",
                        type=int,
                        default=101)
    parser.add_argument('--pos-range', '-pr', help="Number of flanking nucleotides to calculate factors on",
                        type=int,
                        default=12)
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/3_codon_context_score/byAminoAcid_sub/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()

    variant_df = pd.read_csv(args.variant_file,
                             sep="\t",
                             dtype={"CHR":str})
    
    codon_nuc_bias_df = pd.read_csv(args.factor_file,
                                    sep="\t")

    print("Using this scoring table:")
    print(codon_nuc_bias_df.head())
    
    #Add per-base scores
    third_base_pos = int((args.num_pos-1)/2+1)
    variant_scored_df = annotate_variant_context(variant_df,
                                                 codon_nuc_bias_df,
                                                 third_base_pos,
                                                 args.pos_range)
    #Calculate codon context score
    variant_scored_df = calculate_codon_context_score(variant_scored_df,
                                                      args.pos_range)

    
    #For a clean version remove all but summary context score columns
    variant_scored_clean_df = clean_up_score_columns(variant_scored_df,
                                                     args.pos_range)

    #Save data frames
    variant_scored_filename = (args.output_folder+"ssnv_codon_context_score_mi_annotated_"+
                               str(args.pos_range)+"nt"+
                               args.output_tag+".tsv")
    variant_scored_clean_filename = (args.output_folder+"ssnv_codon_context_score_mi_"+
                                     str(args.pos_range)+"nt"+
                                     args.output_tag+".tsv")

    variant_scored_df.to_csv(variant_scored_filename,
                             sep="\t", index=False)
    variant_scored_clean_df.to_csv(variant_scored_clean_filename,
                                   sep="\t", index=False)
