#!/usr/bin/env python

"""

Apply standard normalization of context scores to <grouping column>, 
flag variant sites with gnomAD coverage, and then those observed in gnomAD.
Add additional flag that's randomized within CodonChangeSig
 
Usage: python context_score_normalization_maf.py \
    -vf <scored_variant_file> -s <score column> \
        -g <grouping column> \
            -o <output path>
            
Reads: <scored_variant_file>
Writes: <output path>

"""

import pandas as pd
import numpy as np
import argparse

def apply_zscore_by_group (df, score_column, group_column) :
    """
    Apply standard normaliation to score_column scores in df
    on each subset of group_column
    """
    
    summary_df = df.groupby(group_column).agg({score_column: [np.mean, np.std]})
    print(summary_df.head())
    
    df_xgroup = df.groupby(group_column)
    
    zscore_col = score_column+"_"+group_column+"_zscore"
    df[zscore_col] = 0.0
        
    for group, group_df in df_xgroup :
        score_group_column_zscore = ((group_df[score_column] - 
                                      summary_df.loc[group,(score_column, 'mean')])/
                                     summary_df.loc[group,(score_column, 'std')])
        
        df.loc[group_df.index, zscore_col] = score_group_column_zscore
        
    return df

def gnomad2_3_flag (vdf) :
    """
    Fill in NA values in gnomAD AC and AF columns. Set
    filters for coverage and observed in gnomAD and flag variants 
    that are covered, and covered+observed (y column)
    """

    print("Checking for NA values:")
    print(" 'gnomAD2_EX_AC':", sum(vdf["gnomAD2_EX_AC"].isna()))
    print(" 'gnomAD3_WG_AC':", sum(vdf["gnomAD3_WG_AC"].isna()))

    vdf["gnomAD2_EX_AC"] = vdf["gnomAD2_EX_AC"].fillna(0)
    vdf["gnomAD3_WG_AC"] = vdf["gnomAD3_WG_AC"].fillna(0)
    vdf["gnomAD2_EX_AF"] = vdf["gnomAD2_EX_AF"].fillna(0)
    vdf["gnomAD3_WG_AF"] = vdf["gnomAD3_WG_AF"].fillna(0)

    #Use Coverage and Quality columns set in Notebooks [**]
    gnomad_filter = ((vdf["gnomAD2_Coverage"]) &
                     (vdf["gnomAD3_Coverage"]) &
                     (vdf["gnomAD2_Quality"] == "PASS") &
                     (vdf["gnomAD3_Quality"] == "PASS"))
    gnomad_obs = ((vdf["gnomAD2_EX_AC"] > 0) |
                  (vdf["gnomAD3_WG_AC"] > 0))

    #Set default to "not covered in gnomad"
    vdf["y"] = -1
    #-flag sites covered
    vdf.loc[gnomad_filter, "y"] = 0
    #-flag sites covered and observed
    vdf.loc[(gnomad_filter & gnomad_obs), "y"] = 1
    
    return vdf

def get_CodonChangeSig (r) :
    """
    For variant record, return codon change signature
    """
    
    if r["SNVContext"] == "CpG>TpG" :
        cod_sig = r["REF_Codon"]+"pG>"+r["ALT_Codon"]+"pG"
        #Note: in case of NCG > NCA, the contextualizing CpG is already
        # described by the codon
    else :
        cod_sig = r["REF_Codon"]+">"+r["ALT_Codon"]
        
    return cod_sig

def get_randomized_y (df) :
    """
    randomize y column within each subset of CodonChangeSig
    """
    
    #Add Codon Change Sig
    df["CodonChangeSig"] = df.apply(get_CodonChangeSig,
                                    axis=1)

    df["y_rand"] = -1
    df_xCodonChange = df.groupby("CodonChangeSig")
    for cod_change, cc_df in df_xCodonChange:
        print("For codon change:", cod_change)
        print(cc_df["y"].value_counts(normalize=True))
        
        gnom_in = cc_df["y"] > -1
        gnom_in_index = cc_df[gnom_in].index

        y_rand_values = np.random.permutation(cc_df[gnom_in]["y"].values)
        df.loc[gnom_in_index,"y_rand"] = y_rand_values #make sure this doesn't have indices attached

    return df

if __name__ == "__main__" :

    parser = argparse.ArgumentParser()
    parser.add_argument('--scored-variant-file', '-vf', help="Path to TSV containing sSNV and context score information",
                        required=True)
    parser.add_argument('--score-column', '-s', help="Column name for score to be normalized")
    parser.add_argument('--group-column', '-g', help="Column name for variable to normalize scores over")
    parser.add_argument('--output', '-o', help="Path to store output")
    args = parser.parse_args()
    
    variant_scored_df = pd.read_csv(args.scored_variant_file,
                                    sep="\t",
                                    dtype={"CHR":str})

    #Apply Z-score normalization
    print("Applying Z-score normalization.")
    variant_scored_norm_df = apply_zscore_by_group(variant_scored_df,
                                                   args.score_column,
                                                   args.group_column)
    #Add gnomAD flagging
    print("Flagging SNVs with coverage in gnomAD.")
    variant_scored_norm_gnomad_df = gnomad2_3_flag(variant_scored_norm_df)
    #-add randomized gnomAD-y
    print("Adding a randomized gnomAD observation flag.")
    variant_scored_norm_gnomad_df = get_randomized_y(variant_scored_norm_gnomad_df)
    
    #write output
    variant_scored_norm_df.to_csv(args.output,
                                  sep="\t",
                                  index=False)
