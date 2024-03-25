#!/usr/bin/env python

""" 

Range through percentiles of the score column distribution, for each group in 
group column, fill out contingency table: count # of variants above (below) score threshold 
by # of variants with observed (unobserved) status in gnomAD, and measure 
sensitivity and specificity

Usage: python measure_constraint_thresholds_on_variants_fishers.py \
    -vf <scored variant file> -s <score column> -g <group column> \
        -o <output folder> -t <output tag>
Reads: scored variant file
Writes: <output folder>context_score_threshold_quantile_x<group column>_<output tag>.tsv

"""

import argparse
import pandas as pd
import numpy as np
import scipy.stats

def measure_enrichment_against_thresholds(df,
                                          group_column, 
                                          score_column,
                                          qnum) :
    """
    For each group in group_column, go through qnum quantiles and measure
    sensitivity and specificity w/r/t predicting variants with unobserved status 
    (y = 0)
    """
    
    #Group scores
    df_grouped_ = df.groupby(group_column)
    
    #Output storage
    results_list = []
    for group, group_df in df_grouped_ :
        #Get quantized values
        quantized_scores, bins = pd.qcut(group_df[score_column],
                                   qnum,
                                   labels=False,
                                   retbins=True)
        group_y_counts = group_df["y"].value_counts()
        frac_group = group_y_counts.loc[0]/(group_y_counts.sum())
        #v in bin i => bins[i] < v <= bins[i+1] except for i=0:
        # bins[0] <= v <= bins[i+1]
        
        for i in range(qnum) :
            counts_above = group_df[quantized_scores > i]["y"].value_counts()
            counts_at_below = group_df[quantized_scores <= i]["y"].value_counts()
            
            #Get counts
            if 0 in counts_above :
                counts_above_0 = counts_above.loc[0]
            else :
                counts_above_0 = 0
            if 1 in counts_above :
                counts_above_1 = counts_above.loc[1]
            else :
                counts_above_1 = 0
            if 0 in counts_at_below :
                counts_below_0 = counts_at_below.loc[0]
            else :
                counts_below_0 = 0
            if 1 in counts_at_below :
                counts_below_1 = counts_at_below.loc[1]
            else :
                counts_below_1 = 0
            
            #Count tables
            test_array = np.array([[counts_above_0, counts_above_1],
                                   [counts_below_0, counts_below_1]])
            
            #Test array
            try :
                test_res = scipy.stats.chi2_contingency(test_array,
                                                        correction=False)
            except ValueError :
                print("Test array with 0 entry at percentile",
                      (i+1),
                      ":", test_array)
                test_res = [-1,-1]
            #[0] - test stat
            #[1] - p-value
            #[2] - dof
            #[3] - expected
            
            #Enrichment
            total_above = counts_above.sum()
            total_below = counts_at_below.sum()
            frac0_above = counts_above_0/total_above if total_above > 0 else 0
            frac0_below = counts_below_0/total_below if total_below > 0 else 0
            
            #Sensitivity, specificity
            sensitivity_above = counts_above_0/group_y_counts.loc[0]
            sensitivity_below = counts_below_0/group_y_counts.loc[0]
            specificity_above = 1.-(counts_above_1/group_y_counts.loc[1])
            specificity_below = 1.-(counts_below_1/group_y_counts.loc[1])
            
            result_row_above = {"group":group,
                                "threshold_index":i+1,
                                "threshold_quantile":(i+1)/float(qnum),
                                "threshold_value": bins[i+1],
                                "comparison":"above",
                                "selected_total": total_above,
                                "selected_0":counts_above_0,
                                "selected_1":counts_above_1,
                                "overall_0":group_y_counts.loc[0],
                                "overall_1":group_y_counts.loc[1],
                                "frac0_on_selected":frac0_above,
                                "frac0_on_unselected":frac0_below,
                                "frac0_overall":frac_group,
                                "sensitivity":sensitivity_above,
                                "specificity":specificity_above,
                                "test_stat": test_res[0],
                                "test_pv": test_res[1]}
            result_row_below = {"group":group,
                                "threshold_index":i+1,
                                "threshold_quantile":(i+1)/float(qnum),
                                "threshold_value": bins[i+1],
                                "comparison":"below",
                                "selected_total": total_below,
                                "selected_0":counts_below_0,
                                "selected_1":counts_below_1,
                                "overall_0":group_y_counts.loc[0],
                                "overall_1":group_y_counts.loc[1],
                                "frac0_on_selected":frac0_below,
                                "frac0_on_unselected":frac0_above,
                                "frac0_overall":frac_group,
                                "sensitivity":sensitivity_below,
                                "specificity":specificity_below,
                                "test_stat": test_res[0],
                                "test_pv": test_res[1]}
            
            results_list.append(result_row_above)
            results_list.append(result_row_below)
    #Form table
    results_df = pd.DataFrame(results_list)
    return results_df
    

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--scored-variant-file', '-vf', 
                        help="Path to TSV containing sSNV and Sequence Context information",
                        required=True)
    parser.add_argument('--score-column', '-s', help="Column name for score to be normalized")
    parser.add_argument('--group-column', '-g', type=str,
                        help="Column to group scores by when measuring strength of threshold")
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/3_codon_context_score/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    #Read in scored variant file
    cc_score_df = pd.read_csv(args.scored_variant_file,
                          sep="\t",
                          dtype={"CHR":str},
                          usecols=["CHR", "POS", "REF", "ALT",
                                   "y", args.group_column, args.score_column])
    cc_score_df.head()
        
    #Filter table
    cc_score_gnomad_df = cc_score_df.query("y > -1").copy()
    del cc_score_df
    
    #-threshold set
    number_quantiles = 100
    results_df = measure_enrichment_against_thresholds(cc_score_gnomad_df,
                                                       args.group_column,
                                                       args.score_column,
                                                       number_quantiles)
    
    #Write table
    write_filename = (args.output_folder+
                      "context_score_threshold_quantile_x"+args.group_column+
                      "_"+args.output_tag+".tsv")
    results_df.to_csv(write_filename,
                      sep="\t",
                      index=False)
            
            
    