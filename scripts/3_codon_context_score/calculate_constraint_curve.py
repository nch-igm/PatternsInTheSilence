#!/usr/bin/env python

"""

Split up score column into specified # of quantiles, then for each group
in grouping column, bin scores and assign p(MAF > 0) based on the indicator
column. Write out file to reproduce binned-score  vs. p(MAF > 0) curve for each 
group

Usage: python calculate_constraint_curve.py \
    -vf <scored variant file> -d <score column> \
        -y <indicator columns> -g <grouping column> \
            -b <# of bins> -o <output folder> -t <output tag>
Reads: scored variant file
Writes: <output folder>constraint_curve_<# of bins>bins_<output tag>.tsv

"""

import pandas as pd
import numpy as np
import argparse

def get_binned_score (df, score_column, num_bins) :
    """
    bin score column to specified number of quantiles
    """
    
    df[score_column+"_binned"], score_bins = pd.qcut(df[score_column], 
                                                     num_bins, retbins=True,
                                                     labels=False)
    
    return df, score_bins

def get_constraint_columns (df,
                            y_column,
                            binned_score_column,
                            group_column,
                            score_bins) :
    
    """
    Calculate p(MAF > 0) based on given y_column, over each score bin
    """
    
    #Isolate the columns we need for the MAF plot
    columns = ([binned_score_column] + list(y_column) + 
               ["CHR", group_column])
    vdf_score_columns = df[columns]
    
    agg_dict = {y_col: [np.mean, np.std] for y_col in y_column}
    agg_dict["CHR"] = len
    vdf_score_grouped = vdf_score_columns.\
        groupby([binned_score_column, group_column]).\
        agg(agg_dict)
    
    col_col_names = ([[y_col+"_mean", y_col+"_std"] for y_col in y_column] +
                     ["N"])
    col_names = [x for sublist in col_col_names for x in sublist]
    vdf_score_grouped.columns = col_names
    vdf_score_grouped = vdf_score_grouped.reset_index()
    vdf_score_grouped[binned_score_column+"_left"] = \
        vdf_score_grouped[binned_score_column].apply(lambda x: score_bins[int(x)])
    
    return vdf_score_grouped


if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--scored-variant-file', '-vf', help="Path to TSV containing sSNV and Sequence Context information",
                        required=True)
    parser.add_argument('--score-column', '-s', help="Column containing scores to be binned",
                        default="diff_sum_context_score_REF_Codon_zscore")
    parser.add_argument('--indicator-column', '-y', help="Column(s) containing indicator variable",
                        nargs='*', default="y")
    parser.add_argument('--group-column', '-g', help="Column to group scores along for plotting")
    parser.add_argument('--num-bins', '-b', help="Number of bins to separate score distribution into",
                        type=int, default=100)
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/3_codon_context_score/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    parser.add_argument('--usecols', '-u', help="Select list of columns to import",
                        nargs="+")
    args = parser.parse_args()

    
    if args.usecols is None :
        usecols = ["CHR", "POS", "REF", "ALT", "SNVContext", "CodonChangeSig",
                   "REF_AminoAcid_sub", "ALT_AminoAcid_sub",
                   "REF_Codon", "ALT_Codon",
                   "sum_ref_context_score", "sum_alt_context_score",
                   "diff_sum_context_score", "diff_sum_context_score_REF_Codon_zscore",
                   "y", "y_rand"]
    else :
        usecols = args.usecols
        
    #Read in variant file
    variant_scored_df = pd.read_csv(args.scored_variant_file,
                                    sep="\t",
                                    dtype={"CHR":str},
                                    usecols=usecols)

    variant_scored_gnomad_df = variant_scored_df.query("y > -1").copy()
    
    variant_scored_bin_df, cscore_bins = get_binned_score (variant_scored_gnomad_df, args.score_column, args.num_bins)
    variant_scored_grouped = get_constraint_columns(variant_scored_bin_df, args.indicator_column,
                                                    args.score_column+"_binned", args.group_column,
                                                    cscore_bins)
    
    #Save output columns
    vs_grouped_filename = (args.output_folder+"constraint_curve_"+
                           args.score_column+"_"+
                           str(args.num_bins)+"bins_"+
                           args.output_tag+".tsv")
    
    print("Saving...to", vs_grouped_filename)
    print(variant_scored_grouped.head())
    variant_scored_grouped.to_csv(vs_grouped_filename,
                                  sep="\t",
                                  index=False)
